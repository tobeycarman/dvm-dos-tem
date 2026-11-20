#!/usr/bin/env python
"""Integration tests for stage transitions and restart-file mechanics.

Each successful workflow uses a fresh temporary working directory, the pinned
test-data bundle, and one enabled pixel. The stage chain is:

    pre-run + equilibrium -> spinup -> transient -> scenario

The first stage produces ``restart-eq.nc``; each later stage consumes the
previous stage's restart artifact. The tests verify both the artifacts and the
model's per-pixel ``run_status`` value. Invalid restart paths and invalid
restart stages are tested in separate fresh output directories so expected
failures cannot affect the successful workflow.
"""

import json
import pathlib
import subprocess

import netCDF4 as nc
import pytest

import pyddt.util.runmask
import pyddt.util.setup_working_directory as swd
from pyddt.util.test_data import ensure_bundle


REPO_ROOT = pathlib.Path(__file__).resolve().parents[2]
BUNDLE_LOCK = REPO_ROOT / "testing-data" / "bundle-lock.json"


def _configure_run(run_directory, output_directory, restart_from=""):
    """Set output and restart paths in an isolated run's copied config."""
    config_file = pathlib.Path(run_directory, "config", "config.js")
    with config_file.open(encoding="utf-8") as config_handle:
        config = json.load(config_handle)

    config["IO"]["output_dir"] = str(output_directory)
    config["IO"]["restart_from"] = str(restart_from)
    with config_file.open("w", encoding="utf-8") as config_handle:
        json.dump(config, config_handle, indent=2)
    return config_file


def _run_stage(run_directory, output_directory, arguments, restart_from=""):
    """Run one successful model stage and return its output directory."""
    config_file = _configure_run(
        run_directory, output_directory, restart_from=restart_from
    )
    command = [
        "dvmdostem",
        "-f",
        str(config_file),
        "-l",
        "monitor",
    ] + arguments
    result = subprocess.run(
        command,
        cwd=str(run_directory),
        text=True,
        capture_output=True,
    )
    if result.returncode != 0:
        raise AssertionError(
            "Model stage failed with exit code {}:\n{}".format(
                result.returncode, result.stderr
            )
        )
    return pathlib.Path(output_directory)


def _run_stage_expect_failure(
    run_directory, output_directory, arguments, restart_from
):
    """Run one stage that must reject its restart configuration."""
    config_file = _configure_run(
        run_directory, output_directory, restart_from=restart_from
    )
    command = [
        "dvmdostem",
        "-f",
        str(config_file),
        "-l",
        "monitor",
    ] + arguments
    result = subprocess.run(
        command,
        cwd=str(run_directory),
        text=True,
        capture_output=True,
    )
    assert result.returncode != 0, result.stdout + result.stderr


def _assert_successful_pixel(output_directory):
    """Require the selected pixel to report a completed model stage."""
    status_file = output_directory / "run_status.nc"
    with nc.Dataset(status_file) as dataset:
        assert dataset.variables["run_status"][1, 1] == 100


@pytest.fixture(scope="session")
def bundle_input_directory():
    """Fetch the locked bundle and return its model-input directory."""
    bundle_directory = ensure_bundle(BUNDLE_LOCK)
    input_directory = bundle_directory / "data" / "apex_2km_buffer"
    assert input_directory.is_dir()
    return input_directory


@pytest.fixture
def isolated_run_directory(tmp_path, bundle_input_directory):
    """Create a fresh run workspace and enable only pixel (1, 1)."""
    run_directory = tmp_path / "stage-restart-run"
    swd.cmdline_run(
        swd.cmdline_parse(
            [
                "--force",
                "--input-data-path",
                str(bundle_input_directory),
                str(run_directory),
            ]
        )
    )
    pyddt.util.runmask.cmdline_run(
        pyddt.util.runmask.cmdline_parse(
            [
                "--reset",
                "--yx",
                "1",
                "1",
                str(run_directory / "run-mask.nc"),
            ]
        )
    )
    return run_directory


@pytest.fixture
def stage_outputs(isolated_run_directory):
    """Run and return the complete successful restart chain."""
    run_directory = isolated_run_directory
    pre_run_eq = _run_stage(
        run_directory,
        run_directory / "output_pr_eq",
        ["-p", "5", "-e", "5", "-s", "0", "-t", "0", "-n", "0"],
    )
    spinup = _run_stage(
        run_directory,
        run_directory / "output_sp_from_eq",
        ["-p", "0", "-e", "0", "-s", "5", "-t", "0", "-n", "0"],
        restart_from=pre_run_eq / "restart-eq.nc",
    )
    transient = _run_stage(
        run_directory,
        run_directory / "output_tr_from_sp",
        ["-p", "0", "-e", "0", "-s", "0", "-t", "5", "-n", "0"],
        restart_from=spinup / "restart-sp.nc",
    )
    scenario = _run_stage(
        run_directory,
        run_directory / "output_sc_from_tr",
        ["-p", "0", "-e", "0", "-s", "0", "-t", "0", "-n", "5"],
        restart_from=transient / "restart-tr.nc",
    )
    return {
        "run_directory": run_directory,
        "pr_eq": pre_run_eq,
        "sp": spinup,
        "tr": transient,
        "sc": scenario,
    }


def test_stage_restart_workflow_produces_restart_artifacts(stage_outputs):
    """Every stage writes the restart artifact consumed by the next stage."""
    assert (stage_outputs["pr_eq"] / "restart-eq.nc").exists()
    assert (stage_outputs["sp"] / "restart-sp.nc").exists()
    assert (stage_outputs["tr"] / "restart-tr.nc").exists()
    assert (stage_outputs["sc"] / "restart-sc.nc").exists()


def test_stage_restart_workflow_reports_successful_status(stage_outputs):
    """Every completed stage reports status 100 for the enabled pixel."""
    for stage_name in ("pr_eq", "sp", "tr", "sc"):
        _assert_successful_pixel(stage_outputs[stage_name])


def test_restart_eq_file_cannot_start_pr_or_eq(stage_outputs, tmp_path):
    """A restart file cannot be supplied to the pre-run or equilibrium stage."""
    run_directory = stage_outputs["run_directory"]
    restart_file = stage_outputs["pr_eq"] / "restart-eq.nc"

    _run_stage_expect_failure(
        run_directory,
        tmp_path / "output-invalid-pr-restart",
        ["-p", "5", "-e", "5", "-s", "0", "-t", "0", "-n", "0"],
        restart_from=restart_file,
    )
    _run_stage_expect_failure(
        run_directory,
        tmp_path / "output-invalid-eq-restart",
        ["-p", "0", "-e", "5", "-s", "0", "-t", "0", "-n", "0"],
        restart_from=restart_file,
    )


def test_invalid_restart_path_fails(bundle_input_directory, tmp_path):
    """A missing restart file fails before a stage can run."""
    run_directory = tmp_path / "invalid-restart-run"
    swd.cmdline_run(
        swd.cmdline_parse(
            [
                "--force",
                "--input-data-path",
                str(bundle_input_directory),
                str(run_directory),
            ]
        )
    )
    pyddt.util.runmask.cmdline_run(
        pyddt.util.runmask.cmdline_parse(
            [
                "--reset",
                "--yx",
                "1",
                "1",
                str(run_directory / "run-mask.nc"),
            ]
        )
    )

    missing_restart = run_directory / "missing-restart.nc"
    _run_stage_expect_failure(
        run_directory,
        run_directory / "output-invalid-restart",
        ["-p", "0", "-e", "0", "-s", "5", "-t", "0", "-n", "0"],
        restart_from=missing_restart,
    )
    assert not missing_restart.exists()

#!/bin/bash

#
# bootstrap-sel-custom.sh
#
# Clones SEL repos, sets up dotfiles, (bashrc, vimrc, etc), any other
# user level settings.
#

echo "Running bootstrap-sel-custom.sh script..."
echo "Should not require sudo priviledges to run!"

#
# Note: Not sure if it is better to use tilde or "$HOME" ??
#

#
# add github's key to knownhosts
#
if [[ ! -f ~/.ssh/known_hosts ]]; then
  mkdir -p ~/.ssh
  touch ~/.ssh/known_hosts
fi

echo "Appending github's key to ~/.ssh/known_hosts..."
ssh-keyscan github.com >> "$HOME"/.ssh/known_hosts
chmod 600 "$HOME"/.ssh/known_hosts

#
# Clone our software repos...
#
if [ ! -d "$HOME"/dvm-dos-tem ]
then
  git clone git@github.com:ua-snap/dvm-dos-tem.git "$HOME"/dvm-dos-tem
fi
cd dvm-dos-tem
git remote rename origin upstream
git checkout devel
git pull --ff-only upstream devel:devel
git checkout master
git pull --ff-only upstream master:master
cd ..

#
# setup various preference files (dotfiles)
#

# Not sure how to manage this on Ubuntu automatically, so for now, here are
# the changes I (tbc) have made manually, (8/10/18):
#  - increase HISTSIZE, HISTFILESIZE
#  - uncomment line with "force_color_prompt=yes"
#  - add this at the bottom: 
#    export SITE_SPECIFIC_INCLUDES=-I/usr/include/jsoncpp
#  - enable git branch in prompt according to here:
# https://askubuntu.com/questions/730754/how-do-i-show-the-git-branch-with-colours-in-bash-prompt
# Place this in the .bashrc:
# parse_git_branch() {
#  git branch 2> /dev/null | sed -e '/^[^*]/d' -e 's/* \(.*\)/(\1)/'
# }
# if [ "$color_prompt" = yes ]; then
#  PS1='${debian_chroot:+($debian_chroot)}\[\033[01;32m\]\u@\h\[\033[00m\]:\[\033[01;34m\]\w\[\033[01;31m\] $(parse_git_branch)\[\033[00m\]\$ '
# else
#  PS1='${debian_chroot:+($debian_chroot)}\u@\h:\w$(parse_git_branch)\$ '
# fi
# # THE SIX LINES BELOW are the default prompt and the unset (which were in the original .bashrc)
# #if [ "$color_prompt" = yes ]; then
# #    PS1='${debian_chroot:+($debian_chroot)}\[\033[01;32m\]\u@\h\[\033[00m\]:\[\033[01;34m\]\w\[\033[00m\]\$ '
# #else
# #    PS1='${debian_chroot:+($debian_chroot)}\u@\h:\w\$ '
# #fi
# #unset color_prompt force_color_prompt


# # BASH preferences...
# echo "Setting up bashrc preferences file...."
# cat <<EOF > "$HOME"/.bashrc
# # .bashrc

# # Source global definitions
# if [ -f /etc/bashrc ]; then
#   . /etc/bashrc
# fi

# # Uncomment the following line if you don't like systemctl's auto-paging feature:
# # export SYSTEMD_PAGER=

# # User specific aliases and functions

# # Add git branch to bash prompt...
# source /usr/share/git-core/contrib/completion/git-prompt.sh
# export PS1='[\u@\h \W$(declare -F __git_ps1 &>/dev/null && __git_ps1 " (%s)")]\$ '

# # set up some environment variables
# export SITE_SPECIFIC_INCLUDES=-I/usr/include/jsoncpp
# export PATH=$PATH:/usr/lib64/openmpi/bin
# export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/lib64/openmpi/lib

# EOF

# VIM preferences...
cat <<EOF >> "$HOME"/.vimrc
syntax on              " this is needed to see syntax
set ls=2               " allways show status line
set hlsearch           " highlight searches
set incsearch          " do incremental searching
set ruler              " show the cursor position all the time
set visualbell t_vb=   " turn off error beep/flash
set ignorecase         " ignore case while searching
set number             " put numbers on side
set expandtab          " insert tabs instead of spaces
set tabstop=2          " use 2 spaces
set shiftwidth=2       " how many columns to move with reindent operators (>>, <<)

EOF

# GIT prefs..?
#  - email, editor, color etc??

# Matplotlib prefs...?
#   - might be able to source a file from a gist online?
#   e.g. https://gist.github.com/huyng/816622



echo "DONE setting up a dvm-dos-tem environment. You should be ready to go!"

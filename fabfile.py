#
# Edit `config` line to fit in your environemnt.

# To install fabric and cuisne,
#
#   # update setuptools
#   $ sudo pip install -U setuptools
#   $ sudo pip install setuptools
#
#   $ sudo pip install fabric
#   $ sudo pip install cuisine
#
# You may need to speicfy ARCHFLAGFS on MacOSX environemnt.
# (https://langui.sh/2014/03/10/wunused-command-line-argument-hard-error-in-future-is-a-harsh-mistress/)
#
#   $ sudo ARCHFLAGS=-Wno-error=unused-command-line-argument-hard-error-in-future pip install fabric
#
#
import os, sys
from fabric.api import (env, sudo, put, get, cd, local)
from fabric.utils import puts
from fabric.colors import green
from fabric.decorators import task
from cuisine import (run, dir_ensure, dir_exists, file_exists)

# ----- config --------------------------------------------
env.hosts = ['localhost', 'k', 'linuxbox']
env.use_ssh_config = True

# possible types are 'k_cross', 'linux64' or 'darwin64'
host_type = {
    'k' : 'k_cross'
  , 'xeon.francine' : 'linux64'
  , 'localhost' : 'darwin64'
}

# absolute path to remote build dir.
remote_build_dir = {
    'k' : '/path/to/dir'
  , 'linuxbox' : '/tmp'
  , 'localhost' : '/tmp'
}

# ---------------------------------------------------------

@task
def prepare():
    puts(green('Prepare tools'))
    dir_ensure('deploy')

    if not os.path.exists('./deploy/cmake-2.8.12.2.tar.gz'):
        local('curl http://www.cmake.org/files/v2.8/cmake-2.8.12.2.tar.gz -o deploy/cmake-2.8.12.2.tar.gz')

    put('deploy/cmake-2.8.12.2.tar.gz', remote_build_dir[env.host_string] + '/cmake-2.8.12.2.tar.gz')
    with cd(remote_build_dir[env.host_string]):
        dest_dir = os.path.join(remote_build_dir[env.host_string], 'tools')
        run('tar -zxvf cmake-2.8.12.2.tar.gz')
        run('cd cmake-2.8.12.2; ./configure --prefix=' + dest_dir + ' && make && make install')


@task
def setup():

    if not dir_exists(remote_build_dir[env.host_string]):
        dir_ensure(remote_build_dir[env.host_string], recursive=True)

    dir_ensure(remote_build_dir[env.host_string] + '/build')

    setup_surface()

# Dependency: (None)
@task
def setup_surface():
    puts(green('Configuring SURFACE'))

    local('git archive --format=tar.gz --prefix=SURFACE/ HEAD -o SURFACE.tar.gz')
    put('SURFACE.tar.gz', remote_build_dir[env.host_string] + '/SURFACE.tar.gz')

    with cd(remote_build_dir[env.host_string]):
        run('rm -rf SURFACE')
        run('tar -zxvf SURFACE.tar.gz')
    
    setup_script = ""
    if host_type[env.host_string] == 'k_cross':
        setup_script = './scripts/cmake_k_cross.sh'
    elif host_type[env.host_string] == 'linux64':
        setup_script = './scripts/cmake_linux_x64.sh'
    elif host_type[env.host_string] == 'darwin64':
        setup_script = './scripts/cmake_macosx.sh'
    else:
        print(host_type[env.host_string])
        raise # todo

    with cd(remote_build_dir[env.host_string] + '/SURFACE'):
        run(setup_script)
        run('make -C build')

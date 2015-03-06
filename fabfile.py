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
import json
from fabric.api import (env, sudo, put, get, cd, local)
from fabric.utils import puts
from fabric.colors import green, magenta
from fabric.decorators import task
from cuisine import (run, dir_ensure, dir_exists, file_exists)

# ----- config --------------------------------------------

config_file = "hostconfig.json"
env.use_ssh_config = True

# ---------------------------------------------------------

config = json.loads(open(config_file).read())
env.hosts = config['hosts']

file_hashes = {
    'cmake-2.8.12.2.tar.gz' : '8c6574e9afabcb9fc66f463bb1f2f051958d86c85c37fccf067eb1a44a120e5e'
}


@task
def prepare():
    puts(green('Prepare tools'))
    if not os.path.exists('./deploy'):
        os.mkdir('deploy'):

    remote_build_dir = config['config'][env.host_string]['remote_build_dir']
    build_cmake = config['config'][env.host_string]['build_cmake']
    puts(magenta('remote_build_dir: ' + remote_build_dir))
    puts(magenta('build_cmake: %d' % build_cmake))

    if build_cmake == True:
        if not os.path.exists('./deploy/cmake-2.8.12.2.tar.gz'):
            local('curl http://www.cmake.org/files/v2.8/cmake-2.8.12.2.tar.gz -o deploy/cmake-2.8.12.2.tar.gz')

        put('deploy/cmake-2.8.12.2.tar.gz', remote_build_dir + '/cmake-2.8.12.2.tar.gz')
        with cd(remote_build_dir):
            dest_dir = os.path.join(remote_build_dir, 'tools')
            run('rm -rf cmake-2.8.12.2')
            run('tar -zxvf cmake-2.8.12.2.tar.gz')
            run('cd cmake-2.8.12.2; ./configure --prefix=' + dest_dir + ' && make && make install')


@task
def build():
    remote_build_dir = config['config'][env.host_string]['remote_build_dir']

    if not dir_exists(remote_build_dir):
        dir_ensure(remote_build_dir, recursive=True)

    dir_ensure(remote_build_dir + '/build')

    build_surface()

# Dependency: (None)
@task
def build_surface():
    puts(green('Configuring SURFACE'))
    if not os.path.exists('./deploy'):
        os.mkdir('deploy'):

    remote_build_dir = config['config'][env.host_string]['remote_build_dir']
    host_type = config['config'][env.host_string]['type']
    build_cmake = config['config'][env.host_string]['build_cmake']
    c_compiler = config['config'][env.host_string]['c_compiler']
    cxx_compiler = config['config'][env.host_string]['cxx_compiler']

    local('git archive --format=tar.gz --prefix=SURFACE/ HEAD -o SURFACE.tar.gz')
    put('SURFACE.tar.gz', remote_build_dir + '/SURFACE.tar.gz')

    with cd(remote_build_dir):
        run('rm -rf SURFACE')
        run('tar -zxvf SURFACE.tar.gz')
    
    setup_script = ""
    if host_type == 'k_cross':
        setup_script = './scripts/cmake_k_cross.sh'
    elif host_type == 'linux64':
        setup_script = './scripts/cmake_linux_x64.sh'
    elif host_type == 'darwin64':
        setup_script = './scripts/cmake_macosx.sh'
    else:
        print(host_type)
        raise # todo

    if build_cmake == True:
        cmake_bin_path = os.path.join(config['config'][env.host_string]['remote_build_dir'], "tools/bin/cmake")
        setup_script = "CMAKE_BIN=" + cmake_bin_path + ' ' + setup_script

    if c_compiler is not None:
        setup_script = "CC=" + c_compiler + ' ' + setup_script

    if cxx_compiler is not None:
        setup_script = "CXX=" + cxx_compiler + ' ' + setup_script

    with cd(remote_build_dir + '/SURFACE'):
        run(setup_script)
        run('make -C build')

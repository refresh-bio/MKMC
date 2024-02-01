#!/usr/bin/env python3
import subprocess
import fileinput
import tarfile
import os
import sys

#from https://stackoverflow.com/questions/14697629/running-a-bat-file-though-python-in-current-process
def init_vsvars():
    vswhere_path = r"%ProgramFiles(x86)%/Microsoft Visual Studio/Installer/vswhere.exe"
    vswhere_path = os.path.expandvars(vswhere_path)
    if not os.path.exists(vswhere_path):
        raise EnvironmentError("vswhere.exe not found at: %s", vswhere_path)

    vs_path = os.popen('"{}" -latest -property installationPath'.format(vswhere_path)).read().rstrip()
    vsvars_path = os.path.join(vs_path, "VC\\Auxiliary\\Build\\vcvars64.bat")

    output = os.popen('"{}" && set'.format(vsvars_path)).read()

    for line in output.splitlines():
        pair = line.split("=", 1)
        if(len(pair) >= 2):
            os.environ[pair[0]] = pair[1]


def get_ver():
    with open("mkmc/Version.h") as f:
        for line in f.readlines():
            line = line.strip()
            if "MKMC_VER" in line:
                return line.split("MKMC_VER")[-1].strip().split("\"")[1]
    print("Error: cannot read MKMC_VER")
    sys.exit(1)

def get_os():
    if os.name == 'nt':
        return 'windows'
    elif os.name == 'posix':
        if os.uname()[0] == 'Linux':
            return 'linux'
        elif os.uname()[0] == 'Darwin':
            return 'mac'
        else:
            print("Error: unknown os", os.uname()[0])
            sys.exit(1)
    else:
        print("Error: unknown os.name", os.name)
        sys.exit(1)

def get_hardware():
    if os.name == 'nt':
        return 'x64' # TODO: do a real check and support ARM also...
    elif os.name == 'posix':
        if os.uname()[4] == 'x86_64':
            return 'x64'
        elif os.uname()[4] == 'aarch64' or os.uname()[4] == 'arm64':
            return 'arm64'
        else:
            print("Error: unknown hardware", os.uname()[4])
            sys.exit(1)
    else:
        print("Error: unknown os.name", os.name)
        sys.exit(1)

def run_cmd(cmd):    
    p = subprocess.Popen(cmd, shell=True)
    p.communicate()

system = get_os()
hardware = get_hardware()

ver = get_ver()

print(f"building\n\tVersion: {ver}\n\tOperating system: {system}\n\tHardware: {hardware}")

run_cmd("git submodule init")
run_cmd("git submodule update")

if system == 'windows':
    init_vsvars()
    run_cmd("devenv mkmc.sln /Build \"Release|x64\"")

    with tarfile.open(f"mkmc-{ver}.{system}.{hardware}.tar.gz", "w:gz") as tar:
        #tar.add(source_dir, arcname=os.path.basename(source_dir))
        tar.add("x64\Release\mkmc.exe", arcname="mkmc.exe")
        tar.add("x64\Release\kmc_tools.exe", arcname="kmc_tools.exe")
else:
    run_cmd("make clean")
    run_cmd("make -j")
    run_cmd(f"cd bin; tar -c * | pigz > ../mkmc-{ver}.{system}.{hardware}.tar.gz; cd ..;")

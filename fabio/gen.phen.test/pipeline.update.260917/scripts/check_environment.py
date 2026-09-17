#!/usr/bin/env python3
"""Fail-fast dependency checks and a fingerprint of the installed runtime.

Only local imports/version checks and a temporary PDF/PNG render are performed.
No package is installed, no service is queried, and no cluster job is submitted.
"""
from __future__ import annotations
import argparse
import importlib.metadata
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path
from common import ROOT, digest, read_json, require, sha, write_json

PYTHON_PROBE = '''import importlib, importlib.metadata, json, sys
modules = ['Bio', 'dendropy', 'numpy', 'scipy', 'pkg_resources']
for name in modules: importlib.import_module(name)
from Bio import AlignIO
from scipy.stats import hypergeom, beta
packages = {n: importlib.metadata.version(n) for n in ['biopython','DendroPy','numpy','scipy','setuptools']}
print(json.dumps(dict(executable=sys.executable, version=sys.version.split()[0], packages=packages)))
'''


def executable(value):
    found = shutil.which(str(value))
    require(found is not None, f'Required executable not found: {value}. Create/activate environment.yml first.')
    return str(Path(found).absolute())


def checked_command(command, env=None, cwd=None, timeout=90):
    p = subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True,
                       env=env, cwd=cwd, timeout=timeout)
    require(p.returncode == 0, f'Environment check failed ({command[0]}): {p.stderr.strip() or p.stdout.strip()}')
    return p.stdout.strip(), p.stderr.strip()


def inspect_runtime(python=None, rscript=None, nextflow=None, prefix=None, scheduler=False):
    prefix = Path(prefix).resolve() if prefix else None
    if prefix:
        require((prefix/'conda-meta').is_dir(), f'Not a Conda environment: {prefix}')
    targets = dict(python=python or (prefix/'bin/python' if prefix else sys.executable),
                   Rscript=rscript or (prefix/'bin/Rscript' if prefix else 'Rscript'),
                   nextflow=nextflow or (prefix/'bin/nextflow' if prefix else 'nextflow'),
                   java=prefix/'bin/java' if prefix else 'java')
    binaries = {k:executable(v) for k,v in targets.items()}
    if prefix:
        require(all(Path(p).resolve().is_relative_to(prefix) for p in binaries.values()),
                'Conda mode cannot mix executables from outside its environment')
    env = dict(os.environ, PYTHONDONTWRITEBYTECODE='1', NXF_DISABLE_CHECK_LATEST='true', NXF_OFFLINE='true')
    if prefix:
        env['PATH'] = str(prefix/'bin') + os.pathsep + env.get('PATH','')
        env['PYTHONNOUSERSITE'] = '1'
    stdout, _ = checked_command([binaries['python'],'-c',PYTHON_PROBE],env=env)
    py = json.loads(stdout)
    require(tuple(map(int,py['version'].split('.')[:2])) >= (3,9), 'Python >=3.9 required (recipe uses 3.11)')
    require(int(py['packages']['setuptools'].split('.')[0]) < 82, 'Bundled CAAStools needs pkg_resources/setuptools <82')
    jout, jerr = checked_command([binaries['java'],'-XshowSettings:properties','-version'],env=env)
    jtext = jout + '\n' + jerr
    jhome = re.search(r'^\s*java.home\s*=\s*(.+)$',jtext,re.M)
    jversion = re.search(r'(?:openjdk|java) version "([^"]+)"',jtext)
    require(jhome is not None and jversion is not None, 'Cannot determine Java version/home')
    require(int(jversion.group(1).split('.')[0]) >= 17, 'Java >=17 required')
    if prefix:
        require(Path(jhome.group(1)).resolve().is_relative_to(prefix),'Java home outside the selected Conda environment')
    env['JAVA_HOME'] = env['NXF_JAVA_HOME'] = jhome.group(1)
    nf, nferr = checked_command([binaries['nextflow'],'-version'],env=env)
    nfversion = re.search(r'version\s+([\d.]+)',nf+'\n'+nferr)
    require(nfversion is not None and nfversion.group(1) == '24.04.2',
            'This bundle is tested/pinned to Nextflow 24.04.2; use environment.yml')
    rout,rerr = checked_command([binaries['Rscript'],'--version'],env=env)
    # Verify actual headless graphics support, not only the presence of Rscript.
    with tempfile.TemporaryDirectory(prefix='caas-env-check-') as folder:
        _, warnings = checked_command([binaries['Rscript'],'--vanilla','-e',
                         'pdf("probe.pdf"); plot(1,1); dev.off(); png("probe.png"); plot(1,1); dev.off()'],env=env,cwd=folder)
        require(all((Path(folder)/f).is_file() and (Path(folder)/f).stat().st_size > 0 for f in ('probe.pdf','probe.png')),
                'R PDF/PNG rendering failed: ' + warnings)
    packages = []
    if prefix:
        for path in sorted((prefix/'conda-meta').glob('*.json')):
            meta = read_json(path)
            packages.append({k:meta.get(k) for k in ('name','version','build','subdir','channel','url','sha256','md5')})
    if scheduler:
        for name in ('sbatch','squeue','scancel'): executable(name)
    # No timestamps in the identity: a compatible resume must remain stable.
    identity = dict(binaries=binaries,python=py,java_version=jversion.group(1),java_home=jhome.group(1),
        nextflow_version=nfversion.group(1),R_version=rout or rerr,conda_prefix=str(prefix) if prefix else None,
        installed_conda_packages_sha256=digest(packages) if prefix else None,
        environment_yml_sha256=sha(ROOT/'environment.yml'))
    return dict(status='pass',identity=identity,runtime_sha256=digest(identity),conda_packages=packages,
                checks=['Python imports','pkg_resources','Java','Nextflow','R','headless PDF/PNG'] + (['SLURM commands'] if scheduler else []))


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--conda-prefix',type=Path)
    p.add_argument('--python'); p.add_argument('--rscript'); p.add_argument('--nextflow')
    p.add_argument('--scheduler',action='store_true'); p.add_argument('--output',type=Path)
    a = p.parse_args()
    result = inspect_runtime(a.python,a.rscript,a.nextflow,a.conda_prefix,a.scheduler)
    if a.output: write_json(a.output,result)
    print(json.dumps(dict(status=result['status'],runtime_sha256=result['runtime_sha256'],runtime=result['identity']),indent=2))

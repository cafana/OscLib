# OscLib -- oscillation calculators

## Building Standalone

### Dependencies

#### Required

* `cetmodules`: FNAL CMake toolkit
  + To fetch
```bash
git clone --depth 1 --branch 3.27.03 https://github.com/FNALssi/cetmodules.git
cd cetmodules; mkdir build; cd build;
cmake .. -DCMAKE_INSTALL_PREFIX=$(readlink -f $(uname)); make install;
export cetmodules_ROOT=$(readlink -f $(uname))
```
  + Already fetched: `export cetmodules_ROOT=/path/to/install/prefix`

* ROOT
* Eigen3
* BOOST
* GSL

#### Optional

* STAN

### Build

```bash
git clone git@github.com:cafana/OscLib.git
cd OscLib; mkdir build; cd build
cmake .. -DCMAKE_INSTALL_PREFIX=$(readlink -f $(uname)); make install;
export OscLib_ROOT=$(readlink -f $(uname))
```

* With STAN:

```bash
cmake .. -DOscLib_USE_STAN=ON -DCMAKE_INSTALL_PREFIX=$(readlink -f $(uname)); make install;
```

## Building @ FNAL

### Interactively via ups

* __NOvA Environment__: Update `setup/nova-offline-ups-externals-development`, and `nova-offline-ups-externals-development-prof`. Notify `#cmake`
* __DUNE Environment__: Update `cmake/ups_env_setup.sh`

- `export QUALIFIER=e20:prof` or `c7:debug`, etc
- `export STAN=stan` or `export STAN=stanfree` or `export STAN=stanthread`
- `jenkins/jenkins_build.sh` (or paste parts of it into your terminal)
- set `OSCLIB_LIB` and `OSCLIB_INC` manually to point to what you just built
- rebuild your test release
- Profit!

### NOvA jenkins

#### Build

- Make changes

```
git push
git tag $NEW_TAG_NUMBER
git push --tags
```

- Navigate to https://buildmaster.fnal.gov/buildmaster/view/Nova/job/external/job/osclib_build/ and click "Build Now".
- Wait

#### Deploy

```
wget https://buildmaster.fnal.gov/buildmaster/view/Nova/job/external/job/osclib_collect/lastSuccessfulBuild/artifact/*zip*/archive.zip
unzip archive.zip
mv archive/* .
rm archive.zip

ssh cvmfs${EXPERIMENT}@oasiscfs.fnal.gov
cvmfs_server transaction ${EXPERIMENT}.opensciencegrid.org
# fetch the files you extracted previously to the correct /cvmfs directory
cvmfs_server publish ${EXPERIMENT}.opensciencegrid.org
```

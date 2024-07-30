# OscLib -- oscillation calculators

## How to build and and test your changes locally

- `export QUALIFIER=e20:prof` or `c7:debug`, etc
- `export STAN=stan` or `export STAN=stanfree`
- `jenkins/jenkins_build.sh` (or paste parts of it into your terminal)
- set `OSCLIB_LIB` and `OSCLIB_INC` manually to point to what you just built
- rebuild your test release
- Profit!

## How to build with jenkins

```
git push
git tag $NEW_TAG_NUMBER
git push --tags
```

- Navigate to https://buildmaster.fnal.gov/buildmaster/view/Nova/job/external/job/osclib_build/ and click "Build Now".
- Wait

## How to deploy

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

### For NOvA

Update `setup/nova-offline-ups-externals-development`, and `nova-offline-ups-externals-development-prof`. Notify `#cmake`

### For DUNE

Update `cmake/ups_env_setup.sh`

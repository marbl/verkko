#!/bin/sh
#
#  This will create a release tarball.
#
#  usage: scripts/buildRelease.sh <version>
#
#  It will:
#    1) Fetch submodules (if needed) then check out the correct version.
#    2) Make a directory "verkko-v<version>" and copy all the
#       juicy bits to it, excluding all the putrid .git bits.
#    3) Change Makefile to set VERSION to "verkko release v<version>"
#    4) Create a tar.gz that directory and report the MD5 signature.
#
#  It assumes it is running in the root of a fresh clone of the verkko repo,
#  but will work on existing clones too.  Note that the src/ directory is
#  copied as is, so any leftover crud gets pulled in too (which is why
#  you should build actual releases from a fresh clone).
#
set -e

version=$1

if [ x${version} = x ] ; then
  echo "usage: $0 <version>"
  echo "  where <version> is the decimal version number, e.g., '1.3.1'."
  exit 1
fi
if [ -d verkko-v${version} ] ; then
  echo "ERROR: output directory 'verkko-v${version} exists.  Please remove."
  exit 1
fi
if [ -d verkko-v${version}.tar.gz ] ; then
  echo "ERROR: tarball 'verkko-v${version}.tar.gz exists.  Please remove."
  exit 1
fi


echo "Fetching submodules."

cd src          && git submodule update --init canu           && cd -
cd src/canu/src && git submodule update --init utility        && cd -
cd src/canu/src && git submodule update --init meryl          && cd -
cd src/canu/src && git submodule update --init seqrequester   && cd -
cd src          && git submodule update --init MBG            && cd -
cd src          && git submodule update --init hifioverlapper && cd -
cd src          && git submodule update --init rukki          && cd -

git submodule update --recursive

echo "Building tarball for verkko-v${version}."

mkdir verkko-v${version}

cp -p ./README.licenses verkko-v${version}/
cp -p ./README.md       verkko-v${version}/

for dir in paper src ; do
  rsync -a -f'exclude .git*' -f 'exclude .vscode' $dir verkko-v${version}/
done

cat src/Makefile | sed s/^VERSION.\*/VERSION\ :=\ verkko\ release\ v${version}/ > verkko-v${version}/src/Makefile

tar -cf verkko-v${version}.tar verkko-v${version}
gzip -9v verkko-v${version}.tar

md5sum verkko-v${version}.tar.gz

echo Done.

exit 0

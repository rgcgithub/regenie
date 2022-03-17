#!/usr/bin/env bash
#####
## For use with Docker
help_menu="
Usage regenie_docker.sh OPTIONS

Options:
	--build     create docker image
	--test      test a generated docker image
	--with-bio  compile with Boost Iostreams library
	--with-mkl  compile with MKL library
	--rg-dir    path to the REGENIE source directory
"

# default variables
action=
HAS_BOOST_IOSTREAM=0
MKLROOT=
DFILE=Dockerfile
TEST_SCRIPT=test/test_docker.sh
RGSRC=$(pwd)

while [[ "$#" -gt 0 ]]; do
  case $1 in
    --build) action=build ;;
    --test) action=test ;;
    --with-bio) HAS_BOOST_IOSTREAM=1 ;;
    --with-mkl) MKLROOT=/mkl/ ;;
    --rg-dir) RGSRC="$2"; shift ;;
    -h|--help) action="" ; break ;;
    *) echo "Unknown parameter passed: $1"; echo "$help_menu"; exit 1 ;;
  esac
  shift
done

if [ "$action" = "" ]; then
  echo "$help_menu"; exit 1
fi

if [ ! -f "${RGSRC}/VERSION" ]; then
  echo "must specify REGENIE source directory using '--rg-dir'"; exit 1
fi
cd $RGSRC
RG_VERSION=$(cat VERSION)

if (( HAS_BOOST_IOSTREAM == 1 )); then
  RG_VERSION+=.gz
fi

# create Docker image
if [ "$action" = "build" ]; then
  echo "Building docker image for REGENIE v${RG_VERSION}"

  if (( HAS_BOOST_IOSTREAM == 1 )); then
    echo Compiling with Boost Iostream library
    LIB_BIO=libboost-iostreams-dev
  fi
  if [ "$MKLROOT" != "" ]; then
    echo Compiling with Intel MKL library
    DFILE=Dockerfile_mkl
  fi

	docker build --rm -f ${DFILE} \
		--no-cache --pull \
		--build-arg BOOST_IO=${HAS_BOOST_IOSTREAM} \
		--build-arg LIB_INSTALL=${LIB_BIO} \
		--build-arg LIB_INSTALL2=${LIB_BIO} \
		--build-arg STATIC=1 \
		-t regenie:v${RG_VERSION} .

elif [ "$action" = "test" ]; then
	${TEST_SCRIPT} . "regenie:v${RG_VERSION}" ${HAS_BOOST_IOSTREAM}
fi

cd -

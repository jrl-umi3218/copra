#! /bin/bash
# read modules, build and install them

# force bash to fail if a command fail
set -e

readonly PROGNAME=$(basename $0)
readonly PROGDIR=$(readlink -m $(dirname $0))
readonly ARGS="$@"

readonly CONFIG_FILE_NAME="build_and_install_config"
readonly CONFIG_FILE=$PROGDIR/$CONFIG_FILE_NAME

build_config_dont_exist() {
	[ ! -e $CONFIG_FILE ]
}

create_build_config() {
	cat > $CONFIG_FILE <<EOF
readonly BUILD_DIR="_build"
readonly BUILD_TYPE="Release"
readonly INSTALL_PREFIX="/usr/local"
readonly BUILD_CORE=
readonly BOOST_ROOT_DIR=
EOF
}

load_build_config() {
	source $CONFIG_FILE
}

build_and_install() {
	local path=$1
	(
		mkdir -p $BUILD_DIR/$path
		cd $BUILD_DIR/$path

		cmake $PROGDIR/$path -G"CodeBlocks - Unix Makefiles"\
			-DCMAKE_BUILD_TYPE=$BUILD_TYPE \
			-DCMAKE_INSTALL_PREFIX=$INSTALL_PREFIX \
			-DBOOST_ROOT=$BOOST_ROOT_DIR

		make -j$BUILD_CORE

		if [ -w $INSTALL_PREFIX ]
		then
			make install
		else
			sudo make install
		fi
	)
}

main() {
	if build_config_dont_exist
	then
		create_build_config
		echo "Please modify $CONFIG_FILE and run $PROGNAME again."
		exit 0
	fi

	load_build_config
	for module in 'eigen-qld' 'eigen-quadprog' 'minieigen'
	do
		build_and_install $module
	done
	build_and_install ./
}
main

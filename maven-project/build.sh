#!/bin/sh

die() {

    rc=$?

    echo "ERROR($rc): $@" 
    exit 1

}

cd $(dirname $0)


mvn clean verify || die "mvn build failed!"

OUT_DIR="R-skeleton/inst/java"


rm -rf $OUT_DIR

mkdir -p $OUT_DIR

rc=1

for f in $(find target -name 'extra*.jar'); do
    echo "found: $f"
    cp -v $f $OUT_DIR/ExtraTrees.jar || die "cannot cp $f to $OUR_DIR"
    rc=0
done

ls -l $OUT_DIR

exit $rc


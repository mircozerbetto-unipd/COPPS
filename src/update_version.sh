#!/bin/bash

old_ver=2.1
new_ver=2.2

cat > sed.inp << EOF
s/C++OPPS $old_ver/C++OPPS $new_ver/
s/Version     : $old_ver/Version     : $new_ver/
EOF

mkdir ./tmp

cp -f ./basis.cpp ./tmp/
sed -f ./sed.inp ./tmp/basis.cpp > ./basis.cpp

cp -f ./basis.h ./tmp/
sed -f ./sed.inp ./tmp/basis.h > ./basis.h

cp -f ./constants.h  ./tmp/
sed -f ./sed.inp ./tmp/constants.h > ./constants.h

cp -f ./copps.cpp  ./tmp/
sed -f ./sed.inp ./tmp/copps.cpp > ./copps.cpp

cp ./copps.h  ./tmp/
sed -f ./sed.inp ./tmp/copps.h > ./copps.h

cp -f ./euler.cpp  ./tmp/
sed -f ./sed.inp ./tmp/euler.cpp > ./euler.cpp

cp -f ./euler.h  ./tmp/
sed -f ./sed.inp ./tmp/euler.h > ./euler.h

cp -f ./experimental.cpp ./tmp/
sed -f ./sed.inp ./tmp/experimental.cpp > ./experimental.cpp

cp -f ./experimental.h ./tmp/
sed -f ./sed.inp ./tmp/experimental.h > ./experimental.h

cp -f ./lanczos.cpp  ./tmp/
sed -f ./sed.inp ./tmp/lanczos.cpp > ./lanczos.cpp

cp -f ./lanczos.h  ./tmp/
sed -f ./sed.inp ./tmp/lanczos.h > ./lanczos.h

cp -f ./matrix.cpp ./tmp/
sed -f ./sed.inp ./tmp/matrix.cpp > ./matrix.cpp

cp -f ./matrix.h ./tmp/
sed -f ./sed.inp ./tmp/matrix.h > ./matrix.h

cp -f ./observables.cpp ./tmp/
sed -f ./sed.inp ./tmp/observables.cpp > ./observables.cpp

cp -f ./physics.cpp ./tmp/
sed -f ./sed.inp ./tmp/physics.cpp > ./physics.cpp

cp -f ./physics.h ./tmp/
sed -f ./sed.inp ./tmp/physics.h > ./physics.h

cp ./prep.h ./tmp/
sed -f ./sed.inp ./tmp/prep.h > ./prep.h

cp -f ./relax.cpp ./tmp/
sed -f ./sed.inp ./tmp/relax.cpp > ./relax.cpp

cp -f ./relax.h ./tmp/
sed -f ./sed.inp ./tmp/relax.h > ./relax.h

cp -f ./s3j.cpp ./tmp/
sed -f ./sed.inp ./tmp/s3j.cpp > ./s3j.cpp

cp -f ./s3j.h ./tmp/
sed -f ./sed.inp ./tmp/s3j.h > ./s3j.h

cp -f ./stvec.cpp ./tmp/
sed -f ./sed.inp ./tmp/stvec.cpp > ./stvec.cpp

cp -f ./stvec.h ./tmp/
sed -f ./sed.inp ./tmp/stvec.h > ./stvec.h

cp -f ./tensor4.cpp ./tmp/
sed -f ./sed.inp ./tmp/tensor4.cpp > ./tensor4.cpp

cp -f ./tensor4.h ./tmp/
sed -f ./sed.inp ./tmp/tensor4.h > ./tensor4.h

cp -f ./tensor5.cpp ./tmp/
sed -f ./sed.inp ./tmp/tensor5.cpp > ./tensor5.cpp

cp -f ./tensor5.h ./tmp/
sed -f ./sed.inp ./tmp/tensor5.h > ./tensor5.h

cp -f ./tensor.cpp ./tmp/
sed -f ./sed.inp ./tmp/tensor.cpp > ./tensor.cpp

cp -f ./tensor.h ./tmp/
sed -f ./sed.inp ./tmp/tensor.h > ./tensor.h

cp -f ./types.h ./tmp/
sed -f ./sed.inp ./tmp/types.h > ./types.h

cp -f ./wigner.cpp ./tmp/
sed -f ./sed.inp ./tmp/wigner.cpp > ./wigner.cpp

cp -f ./wigner.h ./tmp/
sed -f ./sed.inp ./tmp/wigner.h > ./wigner.h

rm -rf ./tmp ./sed.inp

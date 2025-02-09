CMDlib=../
here=$PWD

cd $CMDlib
make Z2
cd $here
mv $CMDlib/_Z2.so python/
mv $CMDlib/Z2.py python/Z2backend.py
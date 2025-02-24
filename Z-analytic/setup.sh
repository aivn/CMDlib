CMDlib=../
here=$PWD

cd $CMDlib
make Z2
make Z3
cd $here
mv $CMDlib/_Z2.so $CMDlib/_Z3.so python/
mv $CMDlib/Z2.py python/Z2backend.py
mv $CMDlib/Z3.py python/Z3backend.py
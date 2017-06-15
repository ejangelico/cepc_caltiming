export SYSTEM=slc5
export MOKKA=/workfs/bes/fucd/Mokka/mokka-08-03
export HEPMC=/publicfs/cms/user/huaqiao/convertStdHep/HepMC-2.01.08
export LD_LIBRARY_PATH=$MOKKA/lib/$SYSTEM:$HEPMC/src:$LD_LIBRARY_PATH
export PATH=$MOKKA/bin/$SYSTEM:$PATH


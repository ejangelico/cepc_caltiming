path="$PWD"
cd $path


inputFilePath="../../data/test/"
inputFileName="electron_1GeV_6893"
outputFilePath="../roots/test/"
outputFileName="AnaHit_${inputFileName}.root"

echo \
"
 <marlin>
   
   <execute>
     <processor name=\"MyAnaHit\"/>
     </execute>

     <global>
     <parameter name=\"LCIOInputFiles\">
$(readlink -f ${inputFilePath})/${inputFileName}.slcio
     </parameter>
     <parameter name=\"GearXMLFile\" value=\"$(readlink -f ${inputFilePath})/GearOutput.xml\"/>
     <parameter name=\"MaxRecordNumber\" value=\"1001\"/>
     <parameter name=\"SkipNEvents\" value=\"0\"/>
     <parameter name=\"SupressCheck\" value=\"false\"/>
     <parameter name=\"Verbosity\" options=\"DEBUG0-4,MESSAGE0-4,WARNING0-4,ERROR0-4,SILENT\"> MESSAGE </parameter>
     <parameter name=\"RandomSeed\" value=\"1234567890\" />
     </global>

     <processor name=\"MyAnaHit\" type=\"AnaHit\">
     <parameter name=\"TreeOutputFile\" type=\"StringVec\">
     $(readlink -f ${outputFilePath})/${outputFileName}
     </parameter>
     <parameter name=\"OverwriteFile\" type=\"int\"> 1 </parameter>
     </processor>

     </marlin>
" > tmp_steer/${inputFileName}_reco.xml

echo \
"
cd ..
source env.sh
Marlin $PWD/tmp_steer/${inputFileName}_reco.xml
cd ./steer

" >  tmp_steer/${inputFileName}_reco.sh

chmod +x  tmp_steer/${inputFileName}_reco.sh
#hep_sub  $SimuWorkDir/$OUTPUTDATA/gamma_${energy}GeV_CellSize${CellSize}mm_L200_Si500_W30_${n}_reco.sh -g higgs -o ${energy}GeV_CellSize${CellSize}mm_${n}.out -e ${energy}GeV_CellSize${CellSize}mm_${n}.err
.  tmp_steer/${inputFileName}_reco.sh


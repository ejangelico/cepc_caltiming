path="$PWD"
cd $path

num=3
n=3
cal=26.71
while [ "$n" -le "$num" ]
do

CellSize=5
inum=2
i=2
while [ "$i" -le "$inum" ]
do

case $i in
1) energy=1 ;;
2) energy=5 ;;
3) energy=10 ;;
4) energy=25 ;;
5) energy=50 ;;
6) energy=75 ;;
7) energy=100 ;;
8) energy=125 ;;
9) energy=150 ;;
10) energy=175 ;;
*) echo "i wrong";;
esac

#OUTPUTDATA="H_diphoton"
#OUTPUTDATA="gamma_CellSize"
OUTPUTDATA="gamma"
export SimuWorkDir=$path/simu/
mkdir -p $SimuWorkDir/$OUTPUTDATA/
export tmpDir=$path/tmp_steer/
mkdir -p $tmpDir

echo \
"
 <marlin>
   
   <execute>
     <processor name=\"MyG2CDHGC\"/>
     <processor name=\"MyMarlinArbor\"/>
     <processor name=\"MyBushConnect\"/>
     <processor name=\"MyAnaClu\"/>
     <processor name=\"MyLCIOOutputProcessor\"/>
     </execute>

     <global>
     <parameter name=\"LCIOInputFiles\">
/besfs/groups/higgs/users/zhaoh/cepc/myWorkSpace/SimplifiedGeometry/Generator/data/ScCal/Overlay/gamma_5GeV_CellSize1mm_L30_W28_Sc2/gamma_${energy}GeV_CellSize1mm_L30_W28_Sc2_xchange_${n}.slcio
     </parameter>
     <parameter name=\"GearXMLFile\" value=\"/besfs/groups/higgs/users/zhaoh/cepc/myWorkSpace/SimplifiedGeometry/Generator/data/ScCal/Overlay/gamma_5GeV_CellSize1mm_L30_W28_Sc2/GearOutput.xml\"/>
     <parameter name=\"MaxRecordNumber\" value=\"1001\"/>
     <parameter name=\"SkipNEvents\" value=\"0\"/>
     <parameter name=\"SupressCheck\" value=\"false\"/>
     <parameter name=\"Verbosity\" options=\"DEBUG0-4,MESSAGE0-4,WARNING0-4,ERROR0-4,SILENT\"> MESSAGE </parameter>
     <parameter name=\"RandomSeed\" value=\"1234567890\" />
     </global>

     <processor name=\"MyLCIOOutputProcessor\" type=\"LCIOOutputProcessor\">
     <parameter name=\"LCIOOutputFile\" type=\"string\" >
         gamma_${energy}GeV_CellSize${CellSize}mm_L30_Sc2_x0_reco.slcio
     </parameter>
     <parameter name=\"LCIOWriteMode\" type=\"string\" value=\"WRITE_NEW\"/>
     </processor>

     <processor name=\"MyG2CDHGC\" type=\"G2CDHGC\">
     <parameter name=\"CalibCalo\" type=\"float\"> ${cal} </parameter>
     </processor>

     <processor name=\"MyMarlinArbor\" type=\"MarlinArbor\">
     </processor>

     <processor name=\"MyBushConnect\" type=\"BushConnect\">
     <parameter name=\"FlagDiagnosis\" type=\"int\">0 </parameter>
     <parameter name=\"MCPMIMIC\" type=\"int\">1 </parameter>
     </processor>

	     <processor name=\"MyAnaClu\" type=\"AnaClu\">
	     <parameter name=\"TreeOutputFile\" type=\"StringVec\">
	     AnaClu_gamma_${energy}GeV_CellSize${CellSize}mm_x${n}.root
	     </parameter>
	     <parameter name=\"OverwriteFile\" type=\"int\"> 1 </parameter>
	     </processor>
     </marlin>
" > tmp_steer/gamma_${energy}GeV_CellSize${CellSize}mm_L30_Sc2_x${n}_reco.xml

echo \
"
cd /besfs/groups/higgs/users/zhaoh/cepc/myWorkSpace/ScECAL
source /besfs/groups/higgs/users/zhaoh/cepc/myWorkSpace/ScECAL/env.sh
cd result/tmp
Marlin $PWD/tmp_steer/gamma_${energy}GeV_CellSize${CellSize}mm_L30_Sc2_x${n}_reco.xml
cd /besfs/groups/higgs/users/zhaoh/cepc/myWorkSpace/ScECAL/steer

" >  $SimuWorkDir/$OUTPUTDATA/gamma_${energy}GeV_CellSize${CellSize}mm_L30_Sc2_x${n}_reco.sh

chmod +x  $SimuWorkDir/$OUTPUTDATA/gamma_${energy}GeV_CellSize${CellSize}mm_L30_Sc2_x${n}_reco.sh 
#hep_sub  $SimuWorkDir/$OUTPUTDATA/gamma_${energy}GeV_CellSize${CellSize}mm_L30_Sc2_x${n}_reco.sh -g higgs -o ${energy}GeV_CellSize${CellSize}mm_${n}.out -e ${energy}GeV_CellSize${CellSize}mm_${n}.err
. $SimuWorkDir/$OUTPUTDATA/gamma_${energy}GeV_CellSize${CellSize}mm_L30_Sc2_x${n}_reco.sh 

let "i+=1"
done
let "n+=1"
done
     #<parameter name=\"EcalHitCollections\" type=\"string\"> DigiSiHit </parameter>

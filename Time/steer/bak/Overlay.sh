path="$PWD"
cd $path

num=10
n=10
cal=26.71
while [ "$n" -le "$num" ]
do
o=3
d=$[n-o]

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
	     <processor name=\"MyOverlay\"/>
	     <processor name=\"MyAnaOverlay\"/>
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
	     gamma_${energy}GeV_CellSize${CellSize}mm_L30_Sc2_x${d}_reco.slcio
	     </parameter>
	     <parameter name=\"LCIOWriteMode\" type=\"string\" value=\"WRITE_NEW\"/>
	     </processor>

	     <processor name=\"MyG2CDHGC\" type=\"G2CDHGC\">
	     <parameter name=\"CalibCalo\" type=\"float\"> ${cal} </parameter>
	     </processor>

	     <processor name=\"MyOverlay\" type=\"Overlay\">
	     <parameter name=\"CollectionMap\" type=\"StringVec\">
	     MCParticle	MCParticle
	     SiCalCollection	SiCalCollection
	     </parameter>
	     <parameter name=\"InputFileNames\" type=\"StringVec\">
	     /besfs/groups/higgs/users/zhaoh/cepc/myWorkSpace/ScECAL/result/tmp/gamma_5GeV_CellSize${CellSize}mm_L30_Sc2_x0_reco.slcio
	     </parameter>
	     <parameter name=\"NumberOverlayEvents\" type=\"int\">1 </parameter>
	     </processor>


	     <processor name=\"MyMarlinArbor\" type=\"MarlinArbor\">
	     </processor>

	     <processor name=\"MyBushConnect\" type=\"BushConnect\">
	     <parameter name=\"FlagDiagnosis\" type=\"int\">0 </parameter>
	     <parameter name=\"MCPMIMIC\" type=\"int\">1 </parameter>
	     </processor>

	     <processor name=\"MyAnaOverlay\" type=\"AnaOverlay\">
	     <parameter name=\"TreeOutputFile\" type=\"StringVec\">
	     AnaOverlay_gamma_${energy}GeV_CellSize${CellSize}mm_x${d}.root
	     </parameter>
	     <parameter name=\"OverwriteFile\" type=\"int\"> 1 </parameter>
	     <parameter name=\"Merge\" type=\"int\"> 1 </parameter>
	     </processor>

	     <processor name=\"MyAnaClu\" type=\"AnaClu\">
	     <parameter name=\"TreeOutputFile\" type=\"StringVec\">
	     AnaClu_gamma_${energy}GeV_CellSize${CellSize}mm_x${d}.root
	     </parameter>
	     <parameter name=\"OverwriteFile\" type=\"int\"> 1 </parameter>
	     </processor>

	     </marlin>
	     " > tmp_steer/gamma_${energy}GeV_CellSize${CellSize}mm_L30_Sc2_x${d}_reco.xml

	     echo \
		     "
cd /besfs/groups/higgs/users/zhaoh/cepc/myWorkSpace/ScECAL
source /besfs/groups/higgs/users/zhaoh/cepc/myWorkSpace/ScECAL/env.sh
cd result/tmp
Marlin /besfs/groups/higgs/users/zhaoh/cepc/myWorkSpace/ScECAL/steer/tmp_steer/gamma_${energy}GeV_CellSize${CellSize}mm_L30_Sc2_x${d}_reco.xml
cd /besfs/groups/higgs/users/zhaoh/cepc/myWorkSpace/ScECAL/steer

		     " >  $SimuWorkDir/$OUTPUTDATA/gamma_${energy}GeV_CellSize${CellSize}mm_L30_Sc2_x${d}_reco.sh

		     chmod +x  $SimuWorkDir/$OUTPUTDATA/gamma_${energy}GeV_CellSize${CellSize}mm_L30_Sc2_x${d}_reco.sh 
#hep_sub  $SimuWorkDir/$OUTPUTDATA/gamma_${energy}GeV_CellSize${CellSize}mm_L30_Sc2_x${d}_reco.sh -g higgs -o ${energy}GeV_CellSize${CellSize}mm_${d}.out -e ${energy}GeV_CellSize${CellSize}mm_${d}.err
.  $SimuWorkDir/$OUTPUTDATA/gamma_${energy}GeV_CellSize${CellSize}mm_L30_Sc2_x${d}_reco.sh 

		     let "i+=1"
		     done
		     let "n+=1"
		     done
	     #<parameter name=\"EcalHitCollections\" type=\"string\"> DigiSiHit </parameter>
	     #<parameter name=\"ThresholdsforArborBuilding\" type=\"FloatVec\"> 2 90 50 1.2 </parameter>

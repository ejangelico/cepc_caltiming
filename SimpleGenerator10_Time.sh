#!/bin/bash
path="$PWD"
datapath="${path}/data/gammas/"
cd $path


energy=10 #GeV
randomSeed=$((RANDOM % 10000)) #random 4 digit int
echo ${randomSeed}
output_file="gamma_${energy}GeV_${randomSeed}"
echo ${output_file}
export tmpDir=$path/tmp_steer/
mkdir -p $tmpDir

echo \
"
/generator/generator particleGun
/gun/position 0 0 0 mm
/gun/direction 0. 1. 0.
/gun/energy ${energy} GeV
/gun/momentumSmearing 0.0 GeV
/gun/particle e-
/run/beamOn 1001

exit
" > $tmpDir/event_${output_file}.macro

echo \
"
/Mokka/init/BatchMode true
/Mokka/init/printLevel 0
/Mokka/init/detectorModel CEPC_v1
/Mokka/init/EditGeometry/rmSubDetector all
/Mokka/init/EditGeometry/newSubDetector SiCal

/Mokka/init/dbHost 202.122.37.75
/Mokka/init/user consult
/Mokka/init/dbPasswd consult

/Mokka/init/globalModelParameter world_box_hx 10000
/Mokka/init/globalModelParameter world_box_hy 10000
/Mokka/init/globalModelParameter world_box_hz 25000
/Mokka/init/globalModelParameter SiCalLayerStructure (W:3.,Si:0.5,PCB:2)*200
/Mokka/init/globalModelParameter SiCalZeroThickReset 0

/Mokka/init/globalModelParameter SiCalInnerRadius  1845
/Mokka/init/globalModelParameter SiCalBarrelHalfZ  2245
/Mokka/init/globalModelParameter SiCalEndcapEta1   10000
/Mokka/init/globalModelParameter SiCalEndcapEta2   10000
/Mokka/init/globalModelParameter SiCalBuildBarrel   1
/Mokka/init/globalModelParameter SiCalEndcapOuterR  2500 
/Mokka/init/globalModelParameter SiCalXCellSize  10
/Mokka/init/globalModelParameter SiCalYCellSize  10

/Mokka/init/lcioFilename ${datapath}${output_file}
/Mokka/init/initialMacroFile $tmpDir/event_${output_file}.macro

/Mokka/init/lcioDetailedShowerMode true
/Mokka/init/userInitBool WriteCompleteHepEvt true
/Mokka/init/lcioWriteMode WRITE_NEW
/Mokka/init/lcioStoreCalHitPosition true

#/Mokka/init/dumpG3 false
#/Mokka/init/visumode false
#/Mokka/init/BFactor 1.0

#/Mokka/init/userDeltaOneStep 1e-4 mm
#/Mokka/init/rangeCut 0.1 mm
/Mokka/init/randomSeed ${randomSeed}

#/Mokka/Visu/Detector/ListGeometryTree
#/Mokka/Visu/Detector/DumpGDML SiCal

" > $tmpDir/init_${output_file}.macro

echo \
"
source ./env.sh
Mokka -U $tmpDir/init_${output_file}.macro

" >  $tmpDir/${output_file}.sh

chmod +x  $tmpDir/${output_file}.sh 
.  $tmpDir/${output_file}.sh 
#hep_sub  $tmpDir/${output_file}.sh -g higgs -o ${output_file}.out -e ${output_file}.err
mv GearOutput.xml ${datapath}

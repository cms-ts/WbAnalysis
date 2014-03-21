#!/bin/sh

ulimit -c 0

cd /gpfs/cms/users/schizzi/Wbb2012/CMSSW_5_3_14_patch1
eval `scramv1 runtime -sh`
cd -

cp /gpfs/cms/users/schizzi/Wbb2012/test/demogenanalyzer_cfg.py job.py

echo "fileList = cms.untracked.vstring()" >> job.py
i=1
while [ $i -le $# ]; do
  file=`echo ${!i} | sed -e 's;/gpfs/grid/srm/cms;;'`
  echo "fileList.extend(['"$file"'])" >> job.py
  file=`basename ${!i} | sed -e 's/patTuple/rootTuple/'`
  echo "process.TFileService.fileName = cms.string('"$file"')" >> job.py
  i=$((i+1))
done
echo "process.source.fileNames = fileList" >> job.py

(time cmsRun -j job.xml job.py) > job.log 2>&1

exit

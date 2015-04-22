#!/bin/bash

./submitMC.py --zgamma --o=ZGToLLG --njobs=50 --name=ZGToLLG --inputFolder=/eos/uscms/store/user/bfrancis/SusyNtuples/cms538v1/ZGToLLG_8TeV-madgraph_Summer12_DR53X-PU_RD1_START53_V7N-v1/
./submitMC.py --zgamma --o=WGToLNuG --njobs=50 --name=WGToLNuG --inputFolder=/eos/uscms/store/user/bfrancis/SusyNtuples/cms538v1/WGToLNuG_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_RD1_START53_V7N-v1/
./submitMC.py --zgamma --o=TTGamma --njobs=25 --name=TTGamma --inputFolder=/eos/uscms/store/user/fanxia/SusyNtuples/cms538v1/TTGamma_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_RD1_START53_V7N-v1/

./submitMC.py --zgamma --o=dyJetsToLL --njobs=50 --name=dyJetsToLL --inputFolder=/eos/uscms/store/user/lpcsusystealth/noreplica/ntuples/dyJetsToLL_v1/

#./submitMC.py --zgamma --o=dy1JetsToLL --njobs=50 --name=dy1JetsToLL --inputFolder=/eos/uscms/store/user/lpcsusystealth/noreplica/ntuples/dy1JetsToLL_v1/
#./submitMC.py --zgamma --o=dy2JetsToLL --njobs=50 --name=dy2JetsToLL --inputFolder=/eos/uscms/store/user/lpcsusystealth/noreplica/ntuples/dy2JetsToLL_v1/
#./submitMC.py --zgamma --o=dy3JetsToLL --njobs=50 --name=dy3JetsToLL --inputFolder=/eos/uscms/store/user/lpcsusystealth/noreplica/ntuples/dy3JetsToLL_v1/
#./submitMC.py --zgamma --o=dy4JetsToLL --njobs=50 --name=dy4JetsToLL --inputFolder=/eos/uscms/store/user/lpcsusystealth/noreplica/ntuples/dy4JetsToLL_v1/

#./submitMC.py --zgamma --o=W1JetsToLNu --njobs=50 --name=W1JetsToLNu --inputFolder=/eos/uscms/store/user/lpcsusystealth/noreplica/ntuples/w1JetsToLNu_v1/
#./submitMC.py --zgamma --o=W2JetsToLNu --njobs=50 --name=W2JetsToLNu --inputFolder=/eos/uscms/store/user/lpcsusystealth/noreplica/ntuples/w2JetsToLNu_v1/
./submitMC.py --zgamma --o=W3JetsToLNu --njobs=50 --name=W3JetsToLNu --inputFolder=/eos/uscms/store/user/lpcsusystealth/noreplica/ntuples/w3JetsToLNu_v1/
./submitMC.py --zgamma --o=W4JetsToLNu --njobs=50 --name=W4JetsToLNu --inputFolder=/eos/uscms/store/user/lpcsusystealth/noreplica/ntuples/w4JetsToLNu_v1/

./submitMC.py --zgamma --o=ttJetsHadronic --njobs=50 --name=ttJetsHadronic --inputFolder=/eos/uscms/store/user/lpcsusystealth/noreplica/ntuples/ttJetsHadronic_v1/
./submitMC.py --zgamma --o=ttJetsSemiLep --njobs=25 --name=ttJetsSemiLep --inputFolder=/eos/uscms/store/user/lpcsusystealth/noreplica/ntuples/ttJetsSemiLept_v1/
./submitMC.py --zgamma --o=ttJetsFullLep --njobs=25 --name=ttJetsFullLep --inputFolder=/eos/uscms/store/user/lpcsusystealth/noreplica/ntuples/ttJetsFullLept_v1/

#./submitMC.py --zgamma --o=ttA_2to5 --njobs=25 --name=ttA_2to5 --inputFolder=/eos/uscms/store/user/bfrancis/SusyNtuples/cms538v1/LHE2EDM_WHIZARD_2to5_ttA/

./submitMC.py --zgamma --o=TBar_s --njobs=4 --name=TBar_s --inputFolder=/eos/uscms/store/user/lpcsusystealth/noreplica/ntuples/TBar_s/
./submitMC.py --zgamma --o=TBar_t --njobs=25 --name=TBar_t --inputFolder=/eos/uscms/store/user/lpcsusystealth/noreplica/ntuples/TBar_t/
./submitMC.py --zgamma --o=TBar_tW --njobs=5 --name=TBar_tW --inputFolder=/eos/uscms/store/user/lpcsusystealth/noreplica/ntuples/TBar_tW/
./submitMC.py --zgamma --o=T_s --njobs=4 --name=T_s --inputFolder=/eos/uscms/store/user/lpcsusystealth/noreplica/ntuples/T_s/
./submitMC.py --zgamma --o=T_t --njobs=25 --name=T_t --inputFolder=/eos/uscms/store/user/lpcsusystealth/noreplica/ntuples/T_t/
./submitMC.py --zgamma --o=T_tW --njobs=5 --name=T_tW --inputFolder=/eos/uscms/store/user/lpcsusystealth/noreplica/ntuples/T_tW/

./submitMC.py --zgamma --o=TTWJets --njobs=5 --name=TTWJets --inputFolder=/eos/uscms/store/user/lpcsusystealth/noreplica/ntuples/TTWJets/
./submitMC.py --zgamma --o=TTZJets --njobs=5 --name=TTZJets --inputFolder=/eos/uscms/store/user/lpcsusystealth/noreplica/ntuples/TTZJets/

./submitMC.py --zgamma --o=WW --njobs=100 --name=WW --inputFolder=/eos/uscms/store/user/bfrancis/SusyNtuples/cms538v1/WW_TuneZ2star_8TeV_pythia6_tauola_Summer12_DR53X-PU_S10_START53_V7A-v1/
./submitMC.py --zgamma --o=WZ --njobs=100 --name=WZ --inputFolder=/eos/uscms/store/user/bfrancis/SusyNtuples/cms538v1/WZ_TuneZ2star_8TeV_pythia6_tauola_Summer12_DR53X-PU_S10_START53_V7A-v1/
./submitMC.py --zgamma --o=ZZ --njobs=100 --name=ZZ --inputFolder=/eos/uscms/store/user/3wolf3/SusyNtuples/cms538v1/ZZ_TuneZ2star_8TeV_pythia6_tauola_Summer12_DR53X-PU_S10_START53_V7A-v1/

./submitData.py --zgamma --o=SingleMuA --njobs=25 --inputFolders=/eos/uscms/store/user/lpcsusystealth/noreplica/ntuples/singleMuA_v1/
./submitData.py --zgamma --o=SingleMuB --njobs=25 --inputFolders=/eos/uscms/store/user/lpcsusystealth/noreplica/ntuples/singleMuB_v1/
./submitData.py --zgamma --o=SingleMuC --njobs=25 --inputFolders=/eos/uscms/store/user/lpcsusystealth/noreplica/ntuples/singleMuC/
./submitData.py --zgamma --o=SingleMuD --njobs=25 --inputFolders=/eos/uscms/store/user/lpcsusystealth/noreplica/ntuples/singleMuD_v1/

./submitData.py --zgamma --o=SingleElectronA --njobs=25 --inputFolders=/eos/uscms/store/user/lpcsusystealth/noreplica/ntuples/singleElectronA_v1/
./submitData.py --zgamma --o=SingleElectronB --njobs=25 --inputFolders=/eos/uscms/store/user/lpcsusystealth/noreplica/ntuples/singleElectronB_v1/
./submitData.py --zgamma --o=SingleElectronC --njobs=25 --inputFolders=/eos/uscms/store/user/lpcsusystealth/noreplica/ntuples/singleElectronC_v1/
./submitData.py --zgamma --o=SingleElectronD --njobs=25 --inputFolders=/eos/uscms/store/user/lpcsusystealth/noreplica/ntuples/singleElectronD_v1/

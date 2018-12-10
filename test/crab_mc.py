from WMCore.Configuration import Configuration
config = Configuration()
config.section_('General')
config.General.transferOutputs = True
config.General.transferLogs = True
config.General.requestName = 'processed'
config.General.workArea = 'Fall17'
config.section_('JobType')
config.JobType.psetName = 'HiggsTozz_MiniAOD_mc.py'
config.JobType.pluginName = 'Analysis'
config.JobType.outputFiles = ['roottree_leptons.root']
config.JobType.inputFiles = ['QGL_80X.db']
config.JobType.maxMemoryMB = 2500
#config.JobType.numCores = 8
config.JobType.allowUndistributedCMSSW = True
config.section_('Data')
config.Data.inputDataset = 'datasetpath'
config.Data.outputDatasetTag = 'processed'
config.Data.publishDBS = 'https://cmsdbsprod.cern.ch:8443/cms_dbs_ph_analysis_01_writer/servlet/DBSServlet'
config.Data.publication = False
config.Data.unitsPerJob = 1000
#config.Data.splitting = 'LumiBased'
config.Data.splitting = 'EventAwareLumiBased'
config.Data.outLFNDirBase = '/store/user/reham/MonoHiggs/Fall2017'
config.section_('User')
config.section_('Site')
config.Site.storageSite = 'T2_IT_Bari'

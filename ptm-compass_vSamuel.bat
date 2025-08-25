@echo ** PTM-Compass v1.4 - vSamuel **
@echo ================================

@echo off
REM ==========================================================================
REM Execution of PTM-Compass modules with the new Peak Modelling
REM ==========================================================================

REM Store the path to files passed as argument
set "FILES_PATH=%1"
REM And the name of the fasta
set "FASTA_FILE=%2"

echo Creating directories at %FILES_PATH%...
mkdir "%FILES_PATH%\01-SHIFTSadapter"
mkdir "%FILES_PATH%\02-DuplicateRemover"
mkdir "%FILES_PATH%\03-DMCalibrator"
mkdir "%FILES_PATH%\04-PeakModeller_vSamuel"
mkdir "%FILES_PATH%\05-PeakSelector_vSamuel"
mkdir "%FILES_PATH%\06-PeakInspector"
mkdir "%FILES_PATH%\07-PeakAssignator"
mkdir "%FILES_PATH%\08-PeakFDRer"
mkdir "%FILES_PATH%\09-DM0Solver"
mkdir "%FILES_PATH%\10-ProteinAssigner"
mkdir "%FILES_PATH%\11-TrunkSolver"
mkdir "%FILES_PATH%\12-ProteinAssigner_2"
mkdir "%FILES_PATH%\13-PeakAssignator_2"
mkdir "%FILES_PATH%\14-BinomialSiteListMaker"
mkdir "%FILES_PATH%\15-SiteSolver"
mkdir "%FILES_PATH%\16-ExperimentSeparator"
mkdir "%FILES_PATH%\17-PDMTableMaker"
mkdir "%FILES_PATH%\18-GroupMaker"
mkdir "%FILES_PATH%\19-Joiner"
echo Done!
echo.

echo Starting modules...
echo. 

echo ** 1-SHIFTS adapter **
echo ======================
python SHIFTSadapter.py -i "%FILES_PATH%\data\*.tsv" -o "%FILES_PATH%\01-SHIFTSadapter"
echo.

echo ** 2-Duplicate Remover **
echo =========================
python DuplicateRemover.py -i "%FILES_PATH%\01-SHIFTSadapter\*.feather" -c "%FILES_PATH%\config\params.ini"
move "%FILES_PATH%\01-SHIFTSadapter\*_Unique.feather" "%FILES_PATH%\02-DuplicateRemover"
move "%FILES_PATH%\01-SHIFTSadapter\*SHIFTS_log.txt" "%FILES_PATH%\02-DuplicateRemover"
echo.

echo ** 3-DM Calibrator **
echo =====================
python DMcalibrator.py -i "%FILES_PATH%\02-DuplicateRemover\*_Unique.feather" -c "%FILES_PATH%\config\params.ini"
move "%FILES_PATH%\02-DuplicateRemover\*_calibrated.feather" "%FILES_PATH%\03-DMCalibrator"
move "%FILES_PATH%\02-DuplicateRemover\*Unique_log.txt" "%FILES_PATH%\03-DMCalibrator"
echo.

echo ** 4-Peak Modeller **
echo =====================
python PeakModeller_vSamuel.py -i "%FILES_PATH%\03-DMCalibrator\*_calibrated.feather" -c "%FILES_PATH%\config\params.ini" -t "%FILES_PATH%\config\peak_modelling.tsv"
move "%FILES_PATH%\03-DMCalibrator\*DMTable.feather" "%FILES_PATH%\04-PeakModeller_vSamuel"
move "%FILES_PATH%\03-DMCalibrator\*_DMHistogram.tsv" "%FILES_PATH%\04-PeakModeller_vSamuel"
move "%FILES_PATH%\03-DMCalibrator\PeakModeller_log.txt" "%FILES_PATH%\04-PeakModeller_vSamuel"
echo.

echo ** 5-Peak Selector **
echo =====================
python PeakSelector_vSamuel.py -i "%FILES_PATH%\04-PeakModeller_vSamuel\*.tsv" -c "%FILES_PATH%\config\params.ini" -t "%FILES_PATH%\config\peak_modelling.tsv"
move "%FILES_PATH%\04-PeakModeller_vSamuel\PeakSelector*" "%FILES_PATH%\05-PeakSelector_vSamuel"
move "%FILES_PATH%\04-PeakModeller_vSamuel\*ApexList*" "%FILES_PATH%\05-PeakSelector_vSamuel"
echo.

echo ** 6-Peak Inspector **
echo ======================
echo.

echo ** 7-Peak Assignator **
echo =======================
python PeakAssignator.py -i "%FILES_PATH%\04-PeakModeller_vSamuel\DMTable.feather" -c "%FILES_PATH%\config\params.ini" -a "%FILES_PATH%\05-PeakSelector_vSamuel\Final_ApexList.txt"
move "%FILES_PATH%\04-PeakModeller_vSamuel\*PeakAssignation*" "%FILES_PATH%\07-PeakAssignator"
echo.

echo ** 8-Peak FDR **
echo ================
python PeakFDRer.py -i "%FILES_PATH%\07-PeakAssignator\DMTable_PeakAssignation.feather" -c "%FILES_PATH%\config\params.ini" -a "%FILES_PATH%\05-PeakSelector_vSamuel\Final_ApexList.txt" -e "%FILES_PATH%\config\experimental_table.tsv"
move "%FILES_PATH%\07-PeakAssignator\*FDR*" "%FILES_PATH%\08-PeakFDRer"
echo.

echo ** 9-DM0 Solver **
echo ==================
python DM0Solver.py -i "%FILES_PATH%\08-PeakFDRer\DMTable_PeakAssignation_FDRfiltered.tsv" -c "%FILES_PATH%\config\params.ini" -a "%FILES_PATH%\05-PeakSelector_vSamuel\Final_ApexList.txt"
move "%FILES_PATH%\08-PeakFDRer\*DM0*" "%FILES_PATH%\09-DM0Solver"
echo.

echo ** 10-Protein Assigner **
echo =========================
python ProteinAssigner.py -i "%FILES_PATH%\09-DM0Solver\DMTable_PeakAssignation_FDRfiltered_DM0S.txt" -c "%FILES_PATH%\config\params.ini" -f "%FILES_PATH%\config\%FASTA_FILE%" -o "%FILES_PATH%\10-ProteinAssigner\DMTable_PeakAssignation_FDRfiltered_DM0S_PA.tsv"
echo.

echo ** 11-TrunkSolver **
echo ====================
python TrunkSolver.py -i "%FILES_PATH%\10-ProteinAssigner\DMTable_PeakAssignation_FDRfiltered_DM0S_PA.tsv" -c "%FILES_PATH%\config\params.ini" -f "%FILES_PATH%\config\%FASTA_FILE%"
move "%FILES_PATH%\10-ProteinAssigner\*TrunkSolver*" "%FILES_PATH%\11-TrunkSolver"
move "%FILES_PATH%\10-ProteinAssigner\*_TS.txt" "%FILES_PATH%\11-TrunkSolver"
echo.

echo ** 12-Protein Assigner 2 **
echo ===========================
python ProteinAssigner.py -i "%FILES_PATH%\11-TrunkSolver\DMTable_PeakAssignation_FDRfiltered_DM0S_PA_TS.txt" -c "%FILES_PATH%\config\params_2nd.ini" -f "%FILES_PATH%\config\%FASTA_FILE%" -o "%FILES_PATH%\12-ProteinAssigner_2\DMTable_PeakAssignation_FDRfiltered_DM0S_PA_TS_PA.tsv"
echo.

echo ** 13-Peak Assignator 2 **
echo ==========================
python PeakAssignator.py -i "%FILES_PATH%\12-ProteinAssigner_2\DMTable_PeakAssignation_FDRfiltered_DM0S_PA_TS_PA.tsv" -c "%FILES_PATH%\config\params_2nd.ini" -a "%FILES_PATH%\05-PeakSelector_vSamuel\Final_ApexList.txt"
move "%FILES_PATH%\12-ProteinAssigner_2\DMTable_PeakAssignation_FDRfiltered_DM0S_PA_T_PeakAssignation.tsv" "%FILES_PATH%\13-PeakAssignator_2"
move "%FILES_PATH%\12-ProteinAssigner_2\DMTable_PeakAssignation_FDRfiltered_DM0S_PA_TPeakAssignation_log.txt" "%FILES_PATH%\13-PeakAssignator_2"
echo.

echo ** 14-Binomial Site List **
echo ===========================
python BinomialSiteListMaker.py -i "%FILES_PATH%\13-PeakAssignator_2\DMTable_PeakAssignation_FDRfiltered_DM0S_PA_T_PeakAssignation.tsv" -c "%FILES_PATH%\config\params.ini" -o "%FILES_PATH%\14-BinomialSiteListMaker\BinomialSiteListMaker_PEAKS_Output.xlsx"
echo.

echo ** 15-Site Solver **
echo =====================
python SiteSolver.py -i "%FILES_PATH%\13-PeakAssignator_2\DMTable_PeakAssignation_FDRfiltered_DM0S_PA_T_PeakAssignation.tsv" -c "%FILES_PATH%\config\params.ini" -pl "%FILES_PATH%\config\SiteSolver_List_LabelFree.tsv"
move "%FILES_PATH%\13-PeakAssignator_2\DMTable_PeakAssignation_FDRfiltered_DM0S_PA_T_PeakAssignation_SS.txt" "%FILES_PATH%\15-SiteSolver"
move "%FILES_PATH%\13-PeakAssignator_2\DMTable_PeakAssignation_FDRfiltered_DM0S_PA_T_PeakAssignationSiteSolved_log.txt" "%FILES_PATH%\15-SiteSolver"
echo.

echo ** 16-Experiment Separator **
echo =============================
python ExperimentSeparator.py -i "%FILES_PATH%\15-SiteSolver\DMTable_PeakAssignation_FDRfiltered_DM0S_PA_T_PeakAssignation_SS.txt" -c "Experiment" -o "%FILES_PATH%\16-ExperimentSeparator"
echo.

echo ** 17-PDM Table Maker **
echo ========================
for %%f in ("%FILES_PATH%\16-ExperimentSeparator\*.tsv") do (
	echo %%f > "%FILES_PATH%\16-ExperimentSeparator\experiment.txt"
	python PDMTableMaker.py -i "%FILES_PATH%\16-ExperimentSeparator\experiment.txt" -f "%FILES_PATH%\config\%FASTA_FILE%" -c "%FILES_PATH%\config\params.ini"
	move "%FILES_PATH%\16-ExperimentSeparator\*PDM*" "%FILES_PATH%\17-PDMTableMaker"
	del "%FILES_PATH%\17-PDMTableMaker\experimentPDMTable_log.txt"
)
del "%FILES_PATH%\16-ExperimentSeparator\experiment.txt"
echo.

echo ** 18-GroupMaker **
echo ===================
for %%f in ("%FILES_PATH%\17-PDMTableMaker\*.txt") do (
	python GroupMaker.py -i "%%f" -u "%FILES_PATH%\config\GroupMaker_List_LabelFree.tsv" -c "%FILES_PATH%\config\params.ini"
)
move "%FILES_PATH%\17-PDMTableMaker\*_log.txt" "%FILES_PATH%\18-GroupMaker"
move "%FILES_PATH%\17-PDMTableMaker\*_GM.txt" "%FILES_PATH%\18-GroupMaker"
echo.

echo ** 19-Joiner **
echo ===============
for %%f in ("%FILES_PATH%\18-GroupMaker\*_GM.txt") do (
	python Joiner.py -i "%%f" -c "%FILES_PATH%\config\params.ini"
)
move "%FILES_PATH%\18-GroupMaker\*_Joiner_log.txt" "%FILES_PATH%\19-Joiner"
move "%FILES_PATH%\18-GroupMaker\*_J.txt" "%FILES_PATH%\19-Joiner"
echo.

@pause
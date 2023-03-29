kofam_for_eggNOG=/dss/dssfs02/lwp-dss-0001/u7b03/u7b03-dss-0000/ra72seh/DB/eggNOG_new/kofamscan_for_eggNOG_pred_proteins/
cd ${kofam_for_eggNOG}

for folderNm in `find -type f -print0 | xargs -0 ls | grep 'result_all.txt'`
do
flName=$(echo ${folderNm} | sed 's|/result_all.txt||g' | sed 's|./||g')
fName=${flName}

cat ${folderNm} | sed -e 's/[[:space:]]\{1,\}/\t/' | cut -f1 > ${kofam_for_eggNOG}${fName}_result_tab1; 
cat ${folderNm} | sed -e 's/[[:space:]]\{1,\}/\t/' | cut -f2 | sed -e 's/[[:space:]]\{1,\}/\t/' | cut -f1 >  ${kofam_for_eggNOG}${fName}_result_tab2; 
cat ${folderNm} | sed -e 's/[[:space:]]\{1,\}/\t/' | cut -f2 | sed -e 's/[[:space:]]\{1,\}/\t/' | cut -f2 | sed -e 's/[[:space:]]\{1,\}/\t/' | cut -f1 > ${kofam_for_eggNOG}${fName}_result_tab3; 
cat ${folderNm} | sed -e 's/[[:space:]]\{1,\}/\t/' | cut -f2 | sed -e 's/[[:space:]]\{1,\}/\t/' | cut -f2 | sed -e 's/[[:space:]]\{1,\}/\t/' | cut -f2 | sed -e 's/[[:space:]]\{1,\}/\t/' | cut -f1 >  ${kofam_for_eggNOG}${fName}_result_tab4; 
cat ${folderNm} | sed -e 's/[[:space:]]\{1,\}/\t/' | cut -f2 | sed -e 's/[[:space:]]\{1,\}/\t/' | cut -f2 | sed -e 's/[[:space:]]\{1,\}/\t/' | cut -f2 | sed -e 's/[[:space:]]\{1,\}/\t/' | cut -f2 | sed -e 's/[[:space:]]\{1,\}/\t/' | cut -f1 > ${kofam_for_eggNOG}${fName}_result_tab5; 
cat ${folderNm} | sed -e 's/[[:space:]]\{1,\}/\t/' | cut -f2 | sed -e 's/[[:space:]]\{1,\}/\t/' | cut -f2 | sed -e 's/[[:space:]]\{1,\}/\t/' | cut -f2 | sed -e 's/[[:space:]]\{1,\}/\t/' | cut -f2 | sed -e 's/[[:space:]]\{1,\}/\t/' | cut -f2 | sed -e 's/[[:space:]]\{1,\}/\t/' | cut -f1  > ${kofam_for_eggNOG}${fName}_result_tab6; 
cat ${folderNm} | sed -e 's/[[:space:]]\{1,\}/\t/' | cut -f2 | sed -e 's/[[:space:]]\{1,\}/\t/' | cut -f2 | sed -e 's/[[:space:]]\{1,\}/\t/' | cut -f2 | sed -e 's/[[:space:]]\{1,\}/\t/' | cut -f2 | sed -e 's/[[:space:]]\{1,\}/\t/' | cut -f2 | sed -e 's/[[:space:]]\{1,\}/\t/' | cut -f2  > ${kofam_for_eggNOG}${fName}_result_tab7; 
paste  ${kofam_for_eggNOG}${fName}_result_tab* | tail -n+3 > ${kofam_for_eggNOG}${fName}_result_all_tab_edited; 
rm -f ${kofam_for_eggNOG}${fName}_result_tab*;  

done


kofam_for_Operon=/dss/dssfs02/lwp-dss-0001/u7b03/u7b03-dss-0000/ra72seh/DB/OperonMapper_protein_annotations/kofamscan/
cd ${kofam_for_Operon}

for folderNm in `find -type f -print0 | xargs -0 ls | grep 'result_all.txt'`
do
flName=$(echo ${folderNm} | sed 's|/result_all.txt||g' | sed 's|./||g')
fName=${flName}

cat ${folderNm} | sed -e 's/[[:space:]]\{1,\}/\t/' | cut -f1 > ${kofam_for_Operon}${fName}_result_tab1; 
cat ${folderNm} | sed -e 's/[[:space:]]\{1,\}/\t/' | cut -f2 | sed -e 's/[[:space:]]\{1,\}/\t/' | cut -f1 >  ${kofam_for_Operon}${fName}_result_tab2; 
cat ${folderNm} | sed -e 's/[[:space:]]\{1,\}/\t/' | cut -f2 | sed -e 's/[[:space:]]\{1,\}/\t/' | cut -f2 | sed -e 's/[[:space:]]\{1,\}/\t/' | cut -f1 > ${kofam_for_Operon}${fName}_result_tab3; 
cat ${folderNm} | sed -e 's/[[:space:]]\{1,\}/\t/' | cut -f2 | sed -e 's/[[:space:]]\{1,\}/\t/' | cut -f2 | sed -e 's/[[:space:]]\{1,\}/\t/' | cut -f2 | sed -e 's/[[:space:]]\{1,\}/\t/' | cut -f1 >  ${kofam_for_Operon}${fName}_result_tab4; 
cat ${folderNm} | sed -e 's/[[:space:]]\{1,\}/\t/' | cut -f2 | sed -e 's/[[:space:]]\{1,\}/\t/' | cut -f2 | sed -e 's/[[:space:]]\{1,\}/\t/' | cut -f2 | sed -e 's/[[:space:]]\{1,\}/\t/' | cut -f2 | sed -e 's/[[:space:]]\{1,\}/\t/' | cut -f1 > ${kofam_for_Operon}${fName}_result_tab5; 
cat ${folderNm} | sed -e 's/[[:space:]]\{1,\}/\t/' | cut -f2 | sed -e 's/[[:space:]]\{1,\}/\t/' | cut -f2 | sed -e 's/[[:space:]]\{1,\}/\t/' | cut -f2 | sed -e 's/[[:space:]]\{1,\}/\t/' | cut -f2 | sed -e 's/[[:space:]]\{1,\}/\t/' | cut -f2 | sed -e 's/[[:space:]]\{1,\}/\t/' | cut -f1  > ${kofam_for_Operon}${fName}_result_tab6; 
cat ${folderNm} | sed -e 's/[[:space:]]\{1,\}/\t/' | cut -f2 | sed -e 's/[[:space:]]\{1,\}/\t/' | cut -f2 | sed -e 's/[[:space:]]\{1,\}/\t/' | cut -f2 | sed -e 's/[[:space:]]\{1,\}/\t/' | cut -f2 | sed -e 's/[[:space:]]\{1,\}/\t/' | cut -f2 | sed -e 's/[[:space:]]\{1,\}/\t/' | cut -f2  > ${kofam_for_Operon}${fName}_result_tab7; 
paste  ${kofam_for_Operon}${fName}_result_tab* | tail -n+3 > ${kofam_for_Operon}${fName}_result_all_tab_edited; 
rm -f ${kofam_for_Operon}${fName}_result_tab*;  

done


kofamkoala_proteins=/dss/dssfs02/lwp-dss-0001/u7b03/u7b03-dss-0000/ra72seh/DB/kofamkoala_new_results/
cd ${kofamkoala_proteins}

for kofamkoalaRes in `ls | grep '_proteins_kofam_results.txt$'`
do
flName=$(echo ${kofamkoalaRes} | sed 's|_proteins_kofam_results.txt||g' | sed 's|./||g')
fName=${flName}

cat ${kofamkoalaRes} | sed -e 's/[[:space:]]\{1,\}/\t/' | cut -f1 > ${kofamkoala_proteins}${fName}_result_tab1; 
cat ${kofamkoalaRes} | sed -e 's/[[:space:]]\{1,\}/\t/' | cut -f2 | sed -e 's/[[:space:]]\{1,\}/\t/' | cut -f1 >  ${kofamkoala_proteins}${fName}_result_tab2; 
cat ${kofamkoalaRes} | sed -e 's/[[:space:]]\{1,\}/\t/' | cut -f2 | sed -e 's/[[:space:]]\{1,\}/\t/' | cut -f2 | sed -e 's/[[:space:]]\{1,\}/\t/' | cut -f1 > ${kofamkoala_proteins}${fName}_result_tab3; 
cat ${kofamkoalaRes} | sed -e 's/[[:space:]]\{1,\}/\t/' | cut -f2 | sed -e 's/[[:space:]]\{1,\}/\t/' | cut -f2 | sed -e 's/[[:space:]]\{1,\}/\t/' | cut -f2 | sed -e 's/[[:space:]]\{1,\}/\t/' | cut -f1 >  ${kofamkoala_proteins}${fName}_result_tab4; 
cat ${kofamkoalaRes} | sed -e 's/[[:space:]]\{1,\}/\t/' | cut -f2 | sed -e 's/[[:space:]]\{1,\}/\t/' | cut -f2 | sed -e 's/[[:space:]]\{1,\}/\t/' | cut -f2 | sed -e 's/[[:space:]]\{1,\}/\t/' | cut -f2 | sed -e 's/[[:space:]]\{1,\}/\t/' | cut -f1 > ${kofamkoala_proteins}${fName}_result_tab5; 
cat ${kofamkoalaRes} | sed -e 's/[[:space:]]\{1,\}/\t/' | cut -f2 | sed -e 's/[[:space:]]\{1,\}/\t/' | cut -f2 | sed -e 's/[[:space:]]\{1,\}/\t/' | cut -f2 | sed -e 's/[[:space:]]\{1,\}/\t/' | cut -f2 | sed -e 's/[[:space:]]\{1,\}/\t/' | cut -f2 | sed -e 's/[[:space:]]\{1,\}/\t/' | cut -f1  > ${kofamkoala_proteins}${fName}_result_tab6; 
cat ${kofamkoalaRes} | sed -e 's/[[:space:]]\{1,\}/\t/' | cut -f2 | sed -e 's/[[:space:]]\{1,\}/\t/' | cut -f2 | sed -e 's/[[:space:]]\{1,\}/\t/' | cut -f2 | sed -e 's/[[:space:]]\{1,\}/\t/' | cut -f2 | sed -e 's/[[:space:]]\{1,\}/\t/' | cut -f2 | sed -e 's/[[:space:]]\{1,\}/\t/' | cut -f2  > ${kofamkoala_proteins}${fName}_result_tab7; 
paste  ${kofamkoala_proteins}${fName}_result_tab* | tail -n+3 > ${kofamkoala_proteins}${fName}_result_all_tab_edited; 
rm -f ${kofamkoala_proteins}${fName}_result_tab*;  

done

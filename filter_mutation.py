#post filter
import pandas as pd
import sys

#f = open('E:/bam/funcotated/wgs_sample_list.txt','r')
f = open('C:/Users/User/Downloads/funcotated/TRUP4/trup4.txt','r')
ff = f.readlines()
f.close()
for s in ff:
   s = s.strip()
   s = s+'/'+s
   #header information
   header = pd.read_table('C:/Users/User/Downloads/funcotated/'+s+'.header')
   header = header.values.tolist()
   
   formats =[]
   lists=[]
   for h in header:
      if 'ID=FUNCOTATION' in h[0]:
         funcoCol = h[0].split(':')[1][:-2].split('|')
      elif 'FORMAT=' in h[0]:
         formats.append(h[0].split('=')[2].split(',')[0])
         lists.append([])
   #file = pd.read_table('E:/bam/funcotated/'+s+'.col_noheader.vcf', header=1)
   file = pd.read_csv('E:/bam/funcotated/'+s+'.col_noheader.vcf', header=1,sep='\t',encoding='latin-1')
   file = file[file['FILTER']== 'PASS']
   
   info = file['INFO']
   info = info.str.split(';',expand=True)
   
   AS_FilterStatus=[]
   AS_SB_TABLE=[]
   DP=[]
   ECNT=[]
   FUNCOTATION=[]
   GERMQ=[]
   MBQ=[]
   MFRL=[]
   MMQ=[]
   MPOS=[]
   NALOD=[]
   NLOD=[]
   POPAF=[]
   RPA=[]
   RU=[]
   STR=[]
   STRQ=[]
   TLOD=[]
   info_col = [AS_FilterStatus, AS_SB_TABLE, DP, ECNT, FUNCOTATION, GERMQ,MBQ, MFRL, MMQ, MPOS, NALOD, NLOD, POPAF, RPA, RU,STR, STRQ, TLOD]
   info_column = ['AS_FilterStatus', 'AS_SB_TABLE', 'DP', 'ECNT', 'FUNCOTATION', 'GERMQ','MBQ', 'MFRL', 'MMQ', 'MPOS', 'NALOD', 'NLOD', 'POPAF', 'RPA', 'RU','STR', 'STRQ', 'TLOD']
   nums = list(range(18))
   num=[]
   number=[]

   for i in range(len(info.index)):
      d = info.iloc[i,:].str.split('=')
      for j in range(18):
         if j<=len(d)-1:
            if d[j] is not None:
               if len(d[j]) >=2 :
                  num = info_column.index(d[j][0])
                  info_col[num].append(d[j][1])
                  number.append(num)
                
            
         elif j == 17:
            num = set(nums) - set(number)
            if len(num) !=18:
               for n in list(num):
                  info_col[n].append('No')
               number=[]
   info2 = pd.DataFrame(info_col).T
   info2.columns = info_column
   
   funcotation = info2['FUNCOTATION']
   funcotation = funcotation.str.strip('[')
   funcotation = funcotation.str.strip(']')
   funcotation = funcotation.str.split('|', expand=True)

   funcotation.columns = funcoCol
   ffN = lists
   ffT = lists
   sample =s.split('/')[0]
   if len(file.columns)==11:
      for i in range(len(file.index)):
         f = file.iloc[i,8].split(':')
      
         n=file.iloc[i,9].split(':')
         t = file.iloc[i,10].split(':')
         
         num =[]
         for j in range(len(f)):
            i = formats.index(f[j])
            ffN[i].append(t[j])
            ffT[i].append(t[j])
            num.append(j)
         number = set(list(range(len(lists))))-set(num)
         for k in number:
            ffN[k].append('No')
            ffT[k].append('No')
      ffN = pd.DataFrame(ffN).T
      ffT = pd.DataFrame(ffT).T
      ncol=[]
      tcol=[]
        

      for col in formats:
         ncol.append(sample+'N_'+col)
         tcol.append(sample+'T_'+col)
      ffN.columns = ncol
      ffT.columns = tcol
      
      file = file.drop(['INFO','FORMAT',sample+'T',sample+'N'],axis=True)
      file.reset_index(inplace=True)
      file = file.drop('index',axis=1)
      final = file.join(funcotation)
      info2 = info2.drop('FUNCOTATION',axis=1)
      final = final.join(info2)
      final = final.join(ffN)
      final = final.join(ffT)
      final['POPAF'] =pd.to_numeric(final['POPAF'])
      final[sample+'T_AF'] =pd.to_numeric(final[sample+'T_AF'])
      final[sample+'N_AF'] =pd.to_numeric(final[sample+'N_AF'])
        
      final['gnomAD_exome_AF']=final['gnomAD_exome_AF'].replace([''],0)
      final['gnomAD_exome_AF'].fillna(0)
      af=[]        
      for i in final['gnomAD_exome_AF']:            
         if (str(i).count('_') >=1) :
            val=str(i).split('_%7C_')
            if val[0] == '':
               if srt(i).strip() == '_%7C_':
                  af.append(0)
               else:
                  af.append(float(val[1]))                        
            else:
               af.append(float(val[0]))
         else:
            af.append(float(i))
      final['gnomAD_exome_AF']=af        

      final['gnomAD_genome_AF']=final['gnomAD_genome_AF'].replace([''],0)
      final['gnomAD_genome_AF'].fillna(0)
      af=[]        
      for i in final['gnomAD_genome_AF']:
         print('geome')
         print(i)
         if (str(i).count('_') >=1) :
            val=str(i).split('_%7C_')
            if val[0] == '':
               if str(i).strip() == '_%7C_':
                  af.append(0)
               else:
                  af.append(float(val[1]))
            else:
               af.append(float(val[0]))
         else:
            af.append(float(i))
      final['gnomAD_genome_AF']=af

   else:
      for i in range(len(file.index)):
         f = file.iloc[i,8].split(':')
         t =file.iloc[i,9].split(':')
         
         num=[]
         for j in range(len(f)):
            i=formats.index(f[j])
            ffT[i].append(t[j])
            num.append(i)
         number = set(list(range(len(lists)))) - set(num)
         for k in number:
            ffT[k].append('No')
      ffT = pd.DataFrame(ffT).T
      tcol=[]

      for col in formats:
         tcol.append(sample+'T_'+col)

      ffT.columns = tcol
      file = file.drop(['INFO','FORMAT',sample+'T'], axis=True)
      file.reset_index(inplace=True)
      file = file.drop('index', axis=1)
      final = file.join(funcotation)
      info2 = info2.drop('FUNCOTATION', axis=1)
      final = final.join(info2)
      final = final.join(ffT)
      final['POPAF'] =pd.to_numeric(final['POPAF'])
      final[sample+'T_AF'] =pd.to_numeric(final[sample+'T_AF'])
      final['gnomAD_exome_AF']=final['gnomAD_exome_AF'].replace([''],0)
      final['gnomAD_exome_AF'].fillna(0)
      af=[]        
      for i in final['gnomAD_exome_AF']:
         if (str(i).count('_') >=1) :
            val=str(i).split('_%7C_')
            if val[0] == '':
               if str(i).strip() == '_%7C_':
                  af.append(0)
               else:
                  af.append(float(val[1]))
            else:
               af.append(float(val[0]))
         else:
            af.append(float(i))
      final['gnomAD_exome_AF']=af        

      final['gnomAD_genome_AF']=final['gnomAD_genome_AF'].replace([''],0)
      final['gnomAD_genome_AF'].fillna(0)
      af=[]        
      for i in final['gnomAD_genome_AF']:
         if (str(i).count('_') >=1) :
            val=str(i).split('_%7C_')
                
            if val[0] == '':
               if str(i).strip() == '_%7C_':
                  af.append(0)
               else:
                  af.append(float(val[1]))
            else:
               af.append(float(val[0]))
         else:
            af.append(float(i))
      final['gnomAD_genome_AF']=af
   final = final[final['POPAF'] >= 1.30103] 
   final = final[final['gnomAD_exome_AF'] <= 0.001] 
   final = final[final['gnomAD_genome_AF'] <= 0.005]     
   final.to_excel('E:/bam/funcotated/'+s+'.post_filter.xlsx')
   final['POS']=final['POS'].apply(str)
   final.to_excel('E:/bam/funcotated/'+s+'.post_filter_for_return.xlsx')

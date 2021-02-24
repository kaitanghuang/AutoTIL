####################################################
##
##   This Python-script can be used to calculate the comprehensive scores of TIL subclasses, TIL class and total TIL from "input-geneLengthUnion-v22.csv", "input-geneEnsemblid_Symbol.csv", "input-tilGenes.csv" and the csv file of count of gene expression
##   (Institute) School of Biology and Biological Engineering, South China University of Technology, Guangzhou, China
##   Version 1.0 23.02.2021
##   Needs modules pandas, numpy, os, matplotlib, argparse
##
####################################################

import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import argparse

######命令行参数输入######
# 创建一个 ArgumentParser 对象
parser = argparse.ArgumentParser(description='A comprehensive automatic scoring tool for tumor infiltrating lymphocytes')
parser.add_argument("inputDirectoryPath", type=str, help="enter input directory path that contains all the input files, for example, /home/username/AutoTIL/input")
parser.add_argument("inputCount", type=str, help="enter csv file name of count of gene expression")
parser.add_argument("--tpm", action="store_true", help="output TPM file of gene expression")
parser.add_argument("--log2tpm", action="store_true", help="output log2(TPM+1) file of gene expression")
parser.add_argument("--tilsubclass", action="store_true", help="output comprehensive scores file of 26 TIL subclasses")
parser.add_argument("--til", action="store_true", help="output comprehensive scores file of 4 TIL classes and total TIL")
parser.add_argument("--sample", type=str, help="enter a sample name to be analyzed")
parser.add_argument("--cancer", type=str, help="enter the cancer type of the sample to be analyzed. The optional cancer types include ACC, BLCA, BRCA, CESC, CHOL, COAD, DLBC, ESCA, GBM, HNSC, KICH, KIRC, KIRP, LAML, LGG, LIHC, LUAD, LUSC, MESO, OV, PAAD, PCPG, PRAD, READ, SARC, SKCM, STAD, TGCT, THCA, THYM, UCEC, UCS, UVM")
parser.add_argument('outputDirectoryPath', type=str, help="enter output directory path where all output files will be generated, for example, /home/username/AutoTIL/output")

args = parser.parse_args()

os.chdir(args.outputDirectoryPath)

print("Program starts execution...")
print("Using input directory path: " + args.inputDirectoryPath)
print("Using output directory path: " + args.outputDirectoryPath)

# 1.将输入的基因表达量Count文件（第一列为基因Ensembl ID，其他列为每个样本的基因表达counts数）标准化转换为Transcripts Per Million (TPM)文件
# 读取基因表达count矩阵文件
print("------------------------")
print("Reading count file of gene expression: " + os.path.join(args.inputDirectoryPath, args.inputCount))
df1 = pd.read_csv(os.path.join(args.inputDirectoryPath, args.inputCount), index_col=0)
# 读取基因长度文件
print("Reading gene length file: " + os.path.join(args.inputDirectoryPath, "input-geneLengthUnion-v22.csv"))
lent = pd.read_csv(os.path.join(args.inputDirectoryPath, "input-geneLengthUnion-v22.csv"), index_col=0, squeeze=True)
# 转换为TPM
nldf = df1.div(lent, axis=0)
tpm = (nldf * 1000000) / nldf.sum()

# 2.将标准化后的TPM文件中的基因Ensembl ID转换为基因Symbol，并对同一个基因Symbol对应的不同的基因Ensembl ID表达量取平均值,并输出此时处理后的TPM文件
# 读取基因Ensembl ID对应的基因Symbol数据文件
print("Reading gene ensembl id and symbol annotation file: " + os.path.join(args.inputDirectoryPath, "input-geneEnsemblid_Symbol.csv"))
s = pd.read_csv(os.path.join(args.inputDirectoryPath, "input-geneEnsemblid_Symbol.csv"), squeeze=True, index_col=0)
# index转换
print("Converting count to TPM...")
tpm.index = s[tpm.index]
# 对同一个基因Symbol对应不同的基因Ensembl ID表达量取平均值作为该基因Symbol的表达量
temp = tpm.reset_index()
tpm = temp.groupby(by='Gene_symbol').mean()
# 输出tpm
if args.tpm:
    print('Writing TPM to ' + os.path.join(args.outputDirectoryPath, "res-TPM.csv"))
    tpm.to_csv("res-TPM.csv")
    print("Done!")
    print("------------------------")
else:
    pass

# 3.将上一步处理后的TPM数据转换为Log2(TPM+1)，并输出此时的Log2(TPM+1)文件
print("Converting TPM to log2(TPM+1)...")
tpmLog = tpm + 1
tpmLog = np.log2(tpmLog)
if args.log2tpm:
    print('Writing log2(TPM+1) to ' + os.path.join(args.outputDirectoryPath, "res-log2(TPM+1).csv"))
    tpmLog.to_csv("res-log2(TPM+1).csv")
    print("Done!")
    print("------------------------")
else:
    pass

# 4.利用Log2(TPM+1)，计算26种TIL亚类
# （Act CD4,Act CD8,Tem CD4,Tem CD8,B2M,TAP1,TAP2,HLA-A,HLA-B,HLA-C,HLA-DPA1,HLA-DPB1,
# HLA-E,HLA-F,PD-1,CTLA-4,LAG3,TIGIT,TIM3,PD-L1,PD-L2,CD27,ICOS,IDO1, MDSC,Treg）的综合值，并输出结果文件
# 5.利用26种TIL亚类的综合值，计算4种TIL大类（MHC,CP,EC,SC）的综合值和TIL总的综合值，并输出结果文件
# 提取要计算的样本名
sampleNames = tpmLog.columns.tolist()
# 从制表符分隔的文本文件中读取TIL相关的基因和相应的权重
# 其中class这一列代表基因所属的类：MHC,CP,EC,SC
# 其中subclass这一列代表基因所属的亚类
# gene_name这一列代表具体的基因名
print("Reading TIL genes file: " + os.path.join(args.inputDirectoryPath, "input-tilGenes.csv"))
TILgenesPath = os.path.join(args.inputDirectoryPath, "input-tilGenes.csv")
TILgenesDf = pd.read_csv(TILgenesPath)
temp = TILgenesDf.subclass.tolist()
TILsubclasses = list(set(temp))
TILsubclasses.sort(key=temp.index) # TILsubclasses代表26个基因亚类
# 初始化：MHC(Antigen Processing)抗原呈递大类基因的综合值
# CP(Checkpoints | Immunomodulators)检查点和免疫调节剂大类基因的综合值
# EC(Effector Cells)效应细胞大类基因的综合值
# SC(Suppressor Cells)抑制细胞大类基因的综合值
# AZ(Average Z score)平均Z综合值（是MHC,CP,EC,SC综合值之和）
MHC = []
CP = []
EC = []
SC = []
AZ = []
# 初始化存放所有样本26种基因亚类综合值的dataframe
subclassScore = pd.DataFrame(columns=("sample_id", "subclass", "class", "weight", "subclassScore"))
# 提取表达矩阵中的基因名字
geneNamesExp = tpmLog.index.tolist()
# 提取TIL基因集中的基因名字
geneNamesTIL = TILgenesDf.gene_name.tolist()
# 将TIL基因与表达文件中的基因进行匹配
# 找到表达文件中没有的TIL基因
# 列出那些表达文件中缺失的TIL基因，或者说是与其命名不同的TIL基因，以便研究者修改TIL_genes.txt文件后继续执行该程序
print("Finding TIL genes that don't match log2(TPM+1) genes...")
noMatchGenes = [i for i in geneNamesTIL if i not in geneNamesExp]
# print("表达文件中不匹配的TIL基因有: ", end="")
if len(noMatchGenes) > 0:
    for gene in noMatchGenes:
        print(gene, end=" ")
    print()
    print("The reason is that these genes in log2(TPM+1) are named differently from TIL genes or they are missing in log2(TPM+1).")
    print("It is recommended to execute the program again after modifying the above unmatched genes in the input-tilGenes.csv file to other aliases in the log2(TPM+1).")
else:
    print("There is no TIL gene that does not match log2(TPM+1) genes.")
print("------------------------")
# 计算TIL subclass, TIL class (MHC,CP,EC,SC), 总的TIL的综合值
print("Calculating the comprehensive scores of TIL subclasses...")
for i in range(len(sampleNames)):
  GE = tpmLog.iloc[:, i].values.tolist() #第i个样本的基因表达列
  mGE = np.mean(GE) #均值
  sGE = np.std(GE) #标准差
  temp = tpmLog
  Z1 = (temp.reindex(TILgenesDf.gene_name).iloc[:, i].to_frame() - mGE) / sGE #基因表达矩阵中的TIL基因的Z值
  W1 = TILgenesDf.weight.to_frame() #权重
  MIG = [] #初始化该样本TIL基因亚类的MIG平均Z值list
  WEIGHT = [] #初始化该样本TIL基因亚类的WEIGHT权重list
  
  for gen in TILsubclasses:
    MIG.append(np.nanmean(Z1.loc[(TILgenesDf.subclass==gen).tolist(), :])) #该gen亚类TIL基因的平均Z值
    WEIGHT.append(np.nanmean(W1.loc[(TILgenesDf.subclass==gen).tolist(), :])) #该gen亚类TIL基因的平均权重
    subclassScore = subclassScore.append(pd.DataFrame({"sample_id":[sampleNames[i]], "subclass":[gen], "class":[list(set(TILgenesDf.loc[TILgenesDf.subclass==gen, "class"].values))[0]], "weight":[list(set(TILgenesDf.loc[TILgenesDf.subclass==gen, "weight"].values))[0]],"subclassScore":[np.nanmean(Z1.loc[(TILgenesDf.subclass==gen).tolist(), :])]}), ignore_index=True)

  #一共26个亚类
  WG = [a*b for a,b in zip(MIG, WEIGHT)] #计算亚类Z值与其权重的乘积
  MHC.append(np.mean(WG[0:10])) #计算第i个sample的MHC大类基因的平均的（Z值与其权重的乘积）
  CP.append(np.mean(WG[10:20])) #计算第i个sample的CP大类基因的平均的（Z值与其权重的乘积）
  EC.append(np.mean(WG[20:24])) #计算第i个sample的EC大类基因的平均的（Z值与其权重的乘积）
  SC.append(np.mean(WG[24:26])) #计算第i个sample的SC大类基因的平均的（Z值与其权重的乘积）
  
  AZ.append(MHC[i]+CP[i]+EC[i]+SC[i]) #计算第i个sample的MHC,CP,EC,SC类基因的平均（z值*权重）之和
# 汇总结果并输出
subclassScore["subclassScore_weight"] = subclassScore["subclassScore"] * subclassScore["weight"]
if args.tilsubclass:
    print('Writing the data of the comprehensive scores of TIL subclasses to ' + os.path.join(args.outputDirectoryPath, "res-tilSubclassesScore.csv"))
    subclassScore.to_csv("res-tilSubclassesScore.csv", index=False) # 输出TIL亚类综合值结果文件
    print("Done!")
    print("------------------------")
else:
    pass

print("Calculating the comprehensive scores of TIL classes and total TIL...")
tilScore = pd.DataFrame({"sample_id":sampleNames, "EC_score":EC, "MHC_score":MHC, "CP_score":CP, "SC_score":SC, "TIL_score":AZ})
if args.til:
    print('Writing the data of the comprehensive scores of TIL classes and total TIL to ' + os.path.join(args.outputDirectoryPath, "res-tilScore.csv"))
    tilScore.to_csv("res-tilScore.csv", index=False) # 输出TIL综合值结果文件（包括TIL大类综合值）
    print("Done!")
    print("------------------------")
else:
    pass

# 6.提供要分析的样本ID，可视化评估该样本26种TIL亚类的综合值组成，输出PDF可视化图片
# 利用上一步计算的TIL亚类综合值数据文件
# 筛选特定sample_id的数据
if args.sample:
    print("Analyzing sample " + args.sample + "...")
    sampleId = args.sample
    data = subclassScore.loc[subclassScore.sample_id==sampleId, :]
    temp = data["subclass"].str.cat(data["weight"].astype("str"), sep=" | ").to_frame()
    temp = temp.rename(columns = { "subclass":"subclassName"})
    data = pd.concat([data, temp], axis=1, sort=False)
    data1 = data.loc[data["class"]=="EC", :]
    data2 = data.loc[data["class"]=="MHC", :]
    data3 = data.loc[data["class"]=="CP", :]
    data4 = data.loc[data["class"]=="SC", :]
    data = pd.concat([data1, data2, data3, data4], axis=0, sort=False)
    # 绘图并保存为PDF
    print("Visualizing the comprehensive scores of TIL subclasses...")
    pdf = PdfPages("res-tilSubclassScore" + "(" + sampleId + ").pdf")
    plt.figure(figsize=(25,16),dpi=80)
    plt.bar(data.loc[data["class"]=="EC", "subclassName"], data.loc[data["class"]=="EC", "subclassScore_weight"], label="EC(Effector cells)")
    plt.bar(data.loc[data["class"]=="MHC", "subclassName"], data.loc[data["class"]=="MHC", "subclassScore_weight"], label="MHC(Antigen processing)")
    plt.bar(data.loc[data["class"]=="CP", "subclassName"], data.loc[data["class"]=="CP", "subclassScore_weight"], label="CP(Checkpoints | Immunomodulators)")
    plt.bar(data.loc[data["class"]=="SC", "subclassName"], data.loc[data["class"]=="SC", "subclassScore_weight"], label="SC(Suppressor cells)")
    plt.xticks(data["subclassName"], data["subclassName"], rotation=90)
    plt.legend(fontsize=20)
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)
    plt.title("The comprehensive scores of TIL subclasses of " + sampleId, fontsize=20)
    # print("The figure is being saved as a PDF...")
    pdf.savefig()
    plt.close()
    pdf.close()
    # print("The figure is saved.")
else:
    pass

# 7.最后评估该样本4种TIL大类和TIL总的综合值在TCGA相同癌症大数据中的水平是阳性还是阴性（以中位数为阈值），并输出结果文件
# 利用之前计算的til值tilScore
# 输入要评估的样本ID，及其癌症类型
if args.sample and args.cancer:
    sampleId = args.sample
    print("Comparing with TCGA-" + args.cancer + "...")
    cancer = args.cancer
    # 提取该样本的各个score
    tilScore = tilScore.reset_index().set_index("sample_id")
    score = tilScore.loc[sampleId, :].to_frame()
    score = score.drop("index", axis=0)
    # 存储TCGA33种癌症的EC,MHC,CP,SC,TIL_score的中位数
    medianScore = pd.DataFrame(columns=("cancer", "EC_score", "MHC_score", "CP_score", "SC_score", "TIL_score"))
    medianScore = medianScore.append(pd.DataFrame({"cancer":["ACC"], "EC_score":[0.874764302878283], "MHC_score":[3.2911131023922904], "CP_score":[-0.0360874712838774], "SC_score":[-0.860323708811053], "TIL_score":[3.1388515257643803]}), ignore_index=True)
    medianScore = medianScore.append(pd.DataFrame({"cancer":["BLCA"], "EC_score":[1.1419567181571901], "MHC_score":[3.440459932503895], "CP_score":[-0.2522891595471015], "SC_score":[-1.2503971686042048], "TIL_score":[2.969797220646475]}), ignore_index=True)
    medianScore = medianScore.append(pd.DataFrame({"cancer":["BRCA"], "EC_score":[1.17686627524269], "MHC_score":[3.3165657378337903], "CP_score":[-0.277518830629237], "SC_score":[-1.3597208379607002], "TIL_score":[2.82116374177635]}), ignore_index=True)
    medianScore = medianScore.append(pd.DataFrame({"cancer":["CESC"], "EC_score":[1.29734913233886], "MHC_score":[3.74697766825516], "CP_score":[-0.579520063659792], "SC_score":[-1.37387046261984], "TIL_score":[3.065958180040575]}), ignore_index=True)
    medianScore = medianScore.append(pd.DataFrame({"cancer":["CHOL"], "EC_score":[1.0231199411659149], "MHC_score":[3.510534883115145], "CP_score":[-0.18567700575353], "SC_score":[-1.2064945488976702], "TIL_score":[3.1477842611733298]}), ignore_index=True)
    medianScore = medianScore.append(pd.DataFrame({"cancer":["COAD"], "EC_score":[1.13604091925582], "MHC_score":[3.516110209055325], "CP_score":[-0.2007591913572725], "SC_score":[-1.2367464389114198], "TIL_score":[3.1778936366948702]}), ignore_index=True)
    medianScore = medianScore.append(pd.DataFrame({"cancer":["DLBC"], "EC_score":[1.816181374402485], "MHC_score":[4.38962846754115], "CP_score":[-1.08996422636387], "SC_score":[-2.0166028323986196], "TIL_score":[3.051321364051765]}), ignore_index=True)
    medianScore = medianScore.append(pd.DataFrame({"cancer":["ESCA"], "EC_score":[1.11532873434583], "MHC_score":[3.35178258175646], "CP_score":[-0.30374009214771003], "SC_score":[-1.24691224105493], "TIL_score":[2.84225539188814]}), ignore_index=True)
    medianScore = medianScore.append(pd.DataFrame({"cancer":["GBM"], "EC_score":[1.087584022379835], "MHC_score":[3.27904107352023], "CP_score":[-0.268742935832625], "SC_score":[-1.395233488896995], "TIL_score":[2.6442699595490353]}), ignore_index=True)
    medianScore = medianScore.append(pd.DataFrame({"cancer":["HNSC"], "EC_score":[1.286049065254055], "MHC_score":[3.789659417927045], "CP_score":[-0.516033462704989], "SC_score":[-1.550321978493745], "TIL_score":[2.9881011465537]}), ignore_index=True)
    medianScore = medianScore.append(pd.DataFrame({"cancer":["KICH"], "EC_score":[0.831667006143987], "MHC_score":[3.61622120333633], "CP_score":[-0.14415323332397598], "SC_score":[-0.8386912665137979], "TIL_score":[3.3304104069009304]}), ignore_index=True)
    medianScore = medianScore.append(pd.DataFrame({"cancer":["KIRC"], "EC_score":[1.33623022614065], "MHC_score":[4.00680244052106], "CP_score":[-0.5728969945268225], "SC_score":[-1.654676328579725], "TIL_score":[3.081702604184579]}), ignore_index=True)
    medianScore = medianScore.append(pd.DataFrame({"cancer":["KIRP"], "EC_score":[0.9994506630246995], "MHC_score":[3.58485188003629], "CP_score":[-0.2077038669284385], "SC_score":[-1.2223214877189], "TIL_score":[3.1608683850519945]}), ignore_index=True)
    medianScore = medianScore.append(pd.DataFrame({"cancer":["LAML"], "EC_score":[1.2436940568470902], "MHC_score":[3.0694214286307204], "CP_score":[0.0811840098980151], "SC_score":[-1.62096471539428], "TIL_score":[2.7173707631944]}), ignore_index=True)
    medianScore = medianScore.append(pd.DataFrame({"cancer":["LGG"], "EC_score":[0.9058678005929252], "MHC_score":[2.75296880745821], "CP_score":[-0.0312810233586415], "SC_score":[-0.939684824517714], "TIL_score":[2.6996741447911603]}), ignore_index=True)
    medianScore = medianScore.append(pd.DataFrame({"cancer":["LIHC"], "EC_score":[0.928971836746949], "MHC_score":[3.636676729769], "CP_score":[-0.0598116395202266], "SC_score":[-1.1151178144286398], "TIL_score":[3.28959605788208]}), ignore_index=True)
    medianScore = medianScore.append(pd.DataFrame({"cancer":["LUAD"], "EC_score":[1.33282654402949], "MHC_score":[3.71302914516232], "CP_score":[-0.5576208336019159], "SC_score":[-1.7132560376528898], "TIL_score":[2.74541745130355]}), ignore_index=True)
    medianScore = medianScore.append(pd.DataFrame({"cancer":["LUSC"], "EC_score":[1.30186784647042], "MHC_score":[3.40008059591341], "CP_score":[-0.5257916853886939], "SC_score":[-1.59221931166883], "TIL_score":[2.56739331079583]}), ignore_index=True)
    medianScore = medianScore.append(pd.DataFrame({"cancer":["MESO"], "EC_score":[1.289957566160385], "MHC_score":[3.801200295472275], "CP_score":[-0.4594090337370915], "SC_score":[-1.69358823150823], "TIL_score":[2.8419427945685003]}), ignore_index=True)
    medianScore = medianScore.append(pd.DataFrame({"cancer":["OV"], "EC_score":[1.0595963660396], "MHC_score":[3.259754248388955], "CP_score":[-0.26131220933580956], "SC_score":[-1.240381085535155], "TIL_score":[2.765279462727335]}), ignore_index=True)
    medianScore = medianScore.append(pd.DataFrame({"cancer":["PAAD"], "EC_score":[1.22001383739379], "MHC_score":[3.6585049603217], "CP_score":[-0.31098311902543296], "SC_score":[-1.64333089117375], "TIL_score":[2.94153985266587]}), ignore_index=True)
    medianScore = medianScore.append(pd.DataFrame({"cancer":["PCPG"], "EC_score":[0.9369856332542534], "MHC_score":[3.335351515389365], "CP_score":[-0.11425520972345499], "SC_score":[-1.12922657245737], "TIL_score":[3.049738534489065]}), ignore_index=True)
    medianScore = medianScore.append(pd.DataFrame({"cancer":["PRAD"], "EC_score":[0.962641771380491], "MHC_score":[3.15955754732962], "CP_score":[-0.0789171231999117], "SC_score":[-0.8942750964881021], "TIL_score":[3.12447756485046]}), ignore_index=True)
    medianScore = medianScore.append(pd.DataFrame({"cancer":["READ"], "EC_score":[1.08668126371276], "MHC_score":[3.46690141126374], "CP_score":[-0.160703223583626], "SC_score":[-1.166469278039095], "TIL_score":[3.186315093364795]}), ignore_index=True)
    medianScore = medianScore.append(pd.DataFrame({"cancer":["SARC"], "EC_score":[1.20653127515222], "MHC_score":[3.4236258428435296], "CP_score":[-0.409279822913556], "SC_score":[-1.50892828995891], "TIL_score":[2.70107963231668]}), ignore_index=True)
    medianScore = medianScore.append(pd.DataFrame({"cancer":["SKCM"], "EC_score":[1.06236447802358], "MHC_score":[3.66421775064236], "CP_score":[-0.207809524726952], "SC_score":[-1.29054873524741], "TIL_score":[3.20548756256591]}), ignore_index=True)
    medianScore = medianScore.append(pd.DataFrame({"cancer":["STAD"], "EC_score":[1.22999606373685], "MHC_score":[3.57527597915688], "CP_score":[-0.33395827029685604], "SC_score":[-1.41282172092448], "TIL_score":[3.0519183982837896]}), ignore_index=True)
    medianScore = medianScore.append(pd.DataFrame({"cancer":["TGCT"], "EC_score":[1.3516363203603903], "MHC_score":[3.2235931839289353], "CP_score":[-0.533090960660411], "SC_score":[-1.4330667584822199], "TIL_score":[2.5730821178434855]}), ignore_index=True)
    medianScore = medianScore.append(pd.DataFrame({"cancer":["THCA"], "EC_score":[1.005807941203015], "MHC_score":[3.427084126884215], "CP_score":[-0.1403663061534575], "SC_score":[-1.175972051949305], "TIL_score":[3.1341251344244254]}), ignore_index=True)
    medianScore = medianScore.append(pd.DataFrame({"cancer":["THYM"], "EC_score":[1.51879160890287], "MHC_score":[3.43523774751979], "CP_score":[-0.35296346660095], "SC_score":[-1.42292950580801], "TIL_score":[3.14200008311027]}), ignore_index=True)
    medianScore = medianScore.append(pd.DataFrame({"cancer":["UCEC"], "EC_score":[1.07097392817039], "MHC_score":[3.38760534989593], "CP_score":[-0.336578326219593], "SC_score":[-1.14352987716884], "TIL_score":[2.97568676980492]}), ignore_index=True)
    medianScore = medianScore.append(pd.DataFrame({"cancer":["UCS"], "EC_score":[1.00050218290762], "MHC_score":[2.8481093715120203], "CP_score":[-0.1773276054936505], "SC_score":[-1.08215213579103], "TIL_score":[2.4351195468516353]}), ignore_index=True)
    medianScore = medianScore.append(pd.DataFrame({"cancer":["UVM"], "EC_score":[0.6716399536324211], "MHC_score":[3.3225751093061002], "CP_score":[0.107995699222387], "SC_score":[-0.8128490152653246], "TIL_score":[3.3250887439374455]}), ignore_index=True)
    medianScore = medianScore.reset_index().set_index("cancer")
    medianScore = medianScore.loc[cancer, :].to_frame().drop("index", axis=0) # 筛选特定癌症的各个中位数值数据
    # 评估该样本在TCGA该癌症中的EC,MHC,CP,SC,TIL_score是阳性还是阴性
    # 初始化空dataframe用来存放该样本各个TIL大类的score，并且存放阴性/阳性信息
    npScore = pd.DataFrame(columns=("sample_class", "TIL_class", "score"))
    if float(score.loc["EC_score", :]) > float(medianScore.loc["EC_score", :]):
        npScore = npScore.append(pd.DataFrame({"sample_class":[sampleId], "TIL_class":["EC(+)"], "score":[float(score.loc["EC_score", :])]}))
        npScore = npScore.append(pd.DataFrame({"sample_class":["TCGA-"+cancer], "TIL_class":["EC(+)"], "score":[float(medianScore.loc["EC_score", :])]}))
        # print("EC(+)", end=" ")
    else:
        npScore = npScore.append(pd.DataFrame({"sample_class":[sampleId], "TIL_class":["EC(-)"], "score":[float(score.loc["EC_score", :])]}))
        npScore = npScore.append(pd.DataFrame({"sample_class":["TCGA-"+cancer], "TIL_class":["EC(-)"], "score":[float(medianScore.loc["EC_score", :])]}))
        # print("EC(-)", end=" ")
        
    if float(score.loc["MHC_score", :]) > float(medianScore.loc["MHC_score", :]):
        npScore = npScore.append(pd.DataFrame({"sample_class":[sampleId], "TIL_class":["MHC(+)"], "score":[float(score.loc["MHC_score", :])]}))
        npScore = npScore.append(pd.DataFrame({"sample_class":["TCGA-"+cancer], "TIL_class":["MHC(+)"], "score":[float(medianScore.loc["MHC_score", :])]}))
        # print("MHC(+)", end=" ")
    else:
        npScore = npScore.append(pd.DataFrame({"sample_class":[sampleId], "TIL_class":["MHC(-)"], "score":[float(score.loc["MHC_score", :])]}))
        npScore = npScore.append(pd.DataFrame({"sample_class":["TCGA-"+cancer], "TIL_class":["MHC(-)"], "score":[float(medianScore.loc["MHC_score", :])]}))
        # print("MHC(-)", end=" ")

    if float(score.loc["CP_score", :]) > float(medianScore.loc["CP_score", :]):
        npScore = npScore.append(pd.DataFrame({"sample_class":[sampleId], "TIL_class":["CP(+)"], "score":[float(score.loc["CP_score", :])]}))
        npScore = npScore.append(pd.DataFrame({"sample_class":["TCGA-"+cancer], "TIL_class":["CP(+)"], "score":[float(medianScore.loc["CP_score", :])]}))
        # print("CP(+)", end=" ")
    else:
        npScore = npScore.append(pd.DataFrame({"sample_class":[sampleId], "TIL_class":["CP(-)"], "score":[float(score.loc["CP_score", :])]}))
        npScore = npScore.append(pd.DataFrame({"sample_class":["TCGA-"+cancer], "TIL_class":["CP(-)"], "score":[float(medianScore.loc["CP_score", :])]}))
        # print("CP(-)", end=" ")

    if float(score.loc["SC_score", :]) > float(medianScore.loc["SC_score", :]):
        npScore = npScore.append(pd.DataFrame({"sample_class":[sampleId], "TIL_class":["SC(+)"], "score":[float(score.loc["SC_score", :])]}))
        npScore = npScore.append(pd.DataFrame({"sample_class":["TCGA-"+cancer], "TIL_class":["SC(+)"], "score":[float(medianScore.loc["SC_score", :])]}))
        # print("SC(+)", end=" ")
    else:
        npScore = npScore.append(pd.DataFrame({"sample_class":[sampleId], "TIL_class":["SC(-)"], "score":[float(score.loc["SC_score", :])]}))
        npScore = npScore.append(pd.DataFrame({"sample_class":["TCGA-"+cancer], "TIL_class":["SC(-)"], "score":[float(medianScore.loc["SC_score", :])]}))
        # print("SC(-)", end=" ")

    if float(score.loc["TIL_score", :]) > float(medianScore.loc["TIL_score", :]):
        npScore = npScore.append(pd.DataFrame({"sample_class":[sampleId], "TIL_class":["TIL(+)"], "score":[float(score.loc["TIL_score", :])]}))
        npScore = npScore.append(pd.DataFrame({"sample_class":["TCGA-"+cancer], "TIL_class":["TIL(+)"], "score":[float(medianScore.loc["TIL_score", :])]}))
        # print("TIL(+)", end=" ")
    else:
        npScore = npScore.append(pd.DataFrame({"sample_class":[sampleId], "TIL_class":["TIL(-)"], "score":[float(score.loc["TIL_score", :])]}))
        npScore = npScore.append(pd.DataFrame({"sample_class":["TCGA-"+cancer], "TIL_class":["TIL(-)"], "score":[float(medianScore.loc["TIL_score", :])]}))
        # print("TIL(-)", end=" ")
    # 输出结果文件
    print("Writing the data of TIL classes and total TIL of this sample compared with TCGA-" + args.cancer + " to " + os.path.join(args.outputDirectoryPath, "res-tilScore(" + sampleId + "_and_TCGA-" + cancer + "_median)" + ".csv"))
    npScore.to_csv("res-tilScore(" + sampleId + "_and_TCGA-" + cancer + "_median)" + ".csv", index=False)
    # 提取该样本的阳性/阴性数据
    temp = list(set(npScore["TIL_class"].tolist()))
    temp.sort(key=npScore["TIL_class"].tolist().index)
    # 可视化
    subclasses = temp
    x=np.arange(5)#柱状图在横坐标上的位置
    bar_width=0.45#设置柱状图的宽度
    tick_label = subclasses
    #绘制并列柱状图
    print("Visualizing the data of TIL classes and total TIL of this sample compared with TCGA-" + args.cancer + "...")
    pdf = PdfPages("res-tilScore(" + sampleId + "_and_TCGA-" + cancer + "_median)" + ".pdf")
    plt.figure(figsize=(25,16),dpi=80)
    plt.bar(x, npScore.loc[npScore["sample_class"]==sampleId, "score"], bar_width, color="r", label=sampleId)
    plt.bar(x+bar_width, npScore.loc[npScore["sample_class"]=="TCGA-"+cancer, "score"], bar_width, color="b", label="TCGA-"+cancer)
    plt.xticks(x+bar_width/2,tick_label)
    plt.legend(fontsize=20)
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)
    plt.title("Comparison of comprehensive scores of TIL classes and total TIL of " + sampleId + " with" +  " median of " + "TCGA-" + cancer, fontsize=20)
    # print("The figure is being saved as a PDF...")
    pdf.savefig()
    plt.close()
    pdf.close()
    # print("The figure is saved.")
    print("Done!")
else:
    pass
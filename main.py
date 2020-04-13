from fun import *

# 在运行程序前请先修改fun.py中的全局变量dir和quantity
# 利用pubtator爬取所有数据
scrap()

# 分别将爬取的摘要和注释信息写入文件
writeFile()

# 把摘要和注释信息去重、统一名字，返回pandas数据框
abstract, df, node1 = preprocess()

# 计算所有的term pair并统计他们的权重和IR值，写入文件中
gdPairs = gDMPairs(abstract, df, 'gd')
edgeWeightAndIR(gdPairs, 'gd', node1)

gmPairs = gDMPairs(abstract, df, 'gm')
edgeWeightAndIR(gmPairs, 'gm', node1)

ggPairs = gDMPairs(abstract, df, 'gg')
edgeWeightAndIR(ggPairs, 'gg', node1)

import pandas as pd
from itertools import combinations
import warnings
import spacy
import networkx as nx
import urllib.request
import urllib.error
import time

warnings.filterwarnings('ignore')

dir = 'data'
quantity = '17000'

def scrap():
    f = open("./%s/pmid%s.txt"%(dir, quantity))
    all = open("./%s/all%s.txt"%(dir, quantity), 'a')
    for pmid in f.readlines():
        print("第%d条数据：%s" % (i, pmid.strip()))
        url = r"https://www.ncbi.nlm.nih.gov/research/pubtator-api/publications/export/pubtator?pmids=" \
              + pmid.strip() + r"&concepts=gene,disease,mutation"
        try:
            response = urllib.request.urlopen(url)
            data = response.read().decode('utf-8')
            all.write(data)
            time.sleep(5.8)
        except urllib.error.URLError as e:
            pass
    f.close()
    all.close()

def writeFile():
    gdm = open('./%s/gdm%s.txt' % (dir, quantity), 'w')
    abstract = open('./%s/abstract%s.txt' % (dir, quantity), 'w')
    with open('./%s/all%s.txt' % (dir, quantity)) as f:
        newLines = ''
        absLines = ''
        for line in f.readlines():
            if '|a|' not in line and '|t|' not in line and line is not '\n':
                if line == '\n':
                    line = line.strip()
                id = line.split('\t')[-1].strip()
                if ';' in id:
                    ids = id.split(';')
                    newLine = ''
                    for i in ids:
                        newLine += '\t'.join(line.split('\t')[:-1]+[i]) + '\n'
                    newLines += newLine
                else:
                    newLines += line
            if '|a|' in line:
                absLines += line
        gdm.write(newLines.strip())
        abstract.write(absLines.strip())


def preprocess():
    # abstract数据预处理
    print('./%s/abstract%s.txt'%(dir, quantity))
    abstract = pd.read_csv('./%s/abstract%s.txt'%(dir, quantity),sep='|')
    abstract.columns = ['pmid','a','text']
    abstract.dropna(axis=0, how='any', inplace=True)

    # 去重
    geneDiseaseTable = pd.read_csv("./%s/gdm%s.txt"%(dir, quantity), sep="\t")
    geneDiseaseTable.columns = ["pmid", "start", "end", "name", "genre", "id"]
    df = geneDiseaseTable.drop_duplicates(subset=["pmid", "id"])
    df['name'] = df['name'].str.lower()

    # 去除只出现一次的记录
    ranking = df.loc[:, 'id'].value_counts()
    df.dropna(axis=0, how='any', inplace=True)
    df = df[df['id'].isin(ranking[ranking != 1].index)]

    # 统计分别统计基因、疾病、变异出现的次数
    node1 = df.loc[(df['genre'] == 'Gene'), 'id'].value_counts()
    node2 = df.loc[(df['genre'] == 'DNAMutation') | (df['genre'] == 'ProteinMutation'), 'id'].value_counts()
    node3 = df.loc[(df['genre'] == 'Disease'), 'id'].value_counts()
    for i in node1.index:
        df.loc[df['id'] == i, 'name'] = df[df['id'] == i]['name'].values[0]
    for i in node2.index:
        df.loc[df['id'] == i, 'name'] = df[df['id'] == i]['name'].values[0]
    for i in node3.index:
        df.loc[df['id'] == i, 'name'] = df[df['id'] == i]['name'].values[0]
    name1 = pd.Series([df[df['id'] == i]['name'].values[0] for i in node1.index])
    name2 = pd.Series([df[df['id'] == i]['name'].values[0] for i in node2.index])
    name3 = pd.Series([df[df['id'] == i]['name'].values[0] for i in node3.index])
    node1 = pd.concat([name1, node1.reset_index(drop=True), pd.Series(['Gene']*len(node1))],axis=1,ignore_index=True)
    node2 = pd.concat([name2,node2.reset_index(drop=True), pd.Series(['Mutation']*len(node2))],axis=1,ignore_index=True)
    node3 = pd.concat([name3,node3.reset_index(drop=True), pd.Series(['Disease']*len(node3))],axis=1,ignore_index=True)
    node = pd.concat([node1, node2, node3], axis=0)
    node.to_csv('./%s/gdmNode%s.csv'%(dir, quantity), sep=',', index_label=['index', 'name', 'count', 'genre'])
    return abstract, df, node1


# 计算出现再同一篇摘要中的基因、疾病、变异对作为边
def gDMPairs(abstract, df, type):
    Pairs = []
    pairs = []
    f = open('./%s/sdp%s.txt'%(dir, quantity), 'w')
    for i in df['pmid'].unique():
        gene = df.loc[(df['pmid'] == i) & (df['genre'] == 'Gene') , 'name']
        disease = df.loc[(df['pmid'] == i) & (df['genre'] == 'Disease'), 'name']
        mutation = df.loc[(df['pmid'] == i) & ((df['genre'] == 'ProteinMutation')
                                               | (df['genre'] == 'DNAMutation')), 'name']
        if (abstract.loc[abstract['pmid'] == i, :]).empty:
            continue
        sentences = abstract.loc[abstract['pmid'] == i, 'text'].values[0].split('. ')
        if type == 'gd':
            pairs = [(a, b) for a in gene for b in disease]
        elif type == 'gm':
            pairs = [(a, b) for a in gene for b in mutation]
        elif type == 'gg':
            pairs = list(combinations(gene, 2))
        for pair in pairs:
            x = [(pair[0] in sentence) & (pair[1] in sentence) for sentence in sentences]
            if True in x:
                Pairs.append(list(pair) + [sentences[x.index(True)]])
    f.write('\n'.join(['\t'.join(i) for i in Pairs]))
    return Pairs


# 统计边的权重并输出为csv文件
def edgeWeightAndIR(pairs, kind, node1):
    edgeList = []
    Pairs = [pair[0] for pair in pairs]
    for i in set(Pairs):
        weight = Pairs.count(i)
        edgeList.append(list(i) + [str(weight), str(weight/node1.loc[node1[0] == i[0], 1].values[0]), ])
    csvFile = '\n'.join(['source\ttarget\tweight'] + ['\t'.join(edge) for edge in edgeList])
    with open("./%s/%sNetwork%s.tsv" % (dir, kind, quantity), 'w') as f:
        f.write(csvFile)


nlp = spacy.load("en_core_web_sm")


# 计算pair在一个句子中的shortest dependency path
def SDP(sentence, pair):
    doc = nlp(sentence)
    # for token in doc:
    #     print((token.head.text, token.text, token.dep_))
    edges = []
    for token in doc:
        for child in token.children:
            edges.append(('{0}'.format(token.lower_),
                         '{0}'.format(child.lower_)))
    graph = nx.Graph(edges)
    entity1 = pair[0].lower()
    entity2 = pair[1].lower()
    return nx.shortest_path(graph, source=entity1, target=entity2)


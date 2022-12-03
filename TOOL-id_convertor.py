#!/usr/bin/python
# -*- coding:utf-8 -*- 
########################################
########################################

########################################
### import
import os,re,sys,json,time
import glob
from getopt import getopt,GetoptError
from whoosh.index import open_dir
from whoosh.fields import *
from whoosh.analysis import StemmingAnalyzer
from whoosh.qparser import QueryParser
from whoosh.query import FuzzyTerm,Regex
from whoosh import scoring,qparser
#from whoosh import qparser
########################################
########################################
## adjust the logical for gene search
########################################
########################################

def usage():
    print(
"""
A flexible demo tool for searching Gene Id, Gene Symbol, Ensembl Id , Uniport Id , Refseq Id of protein and Transcript of Gencode.v31(GRCh38) and Ensembl.75(GRCh37).


Usage: python3 path/{0} -i <Search id> [-c Symbol,Uniport_id] [-a -o <FILE>]

Option:
    -a/--all                    Write|Print all ids to output. 
    -p/--print                  Print some results to the screan when input is a file.
    -v/--version                Version.
    -h/--help                   Print this help menu.

Parameter:
    -i/--input  <INT|STR|FILE>  Search id or file, the file should be separated by "\\t" or "," . Default: [150].
                                File will be processed by column.
                                The "xlsx" format is not supported for this version, and will be supported next.

    -c/--custom <STR>           Write|Print the specified id(s) (Symbol,Uniport_id,...). Choose one or more from {1},
                                {2}
                                {11}.
                                Default: Symbol.  
                                "{12}" stands for Transcript id from NCBI(GRCh38). 
                                "{13}" stands for Protein Refseq id from UniProt(GRCh38).
                                "{14}" stands for Uniport id from The Universal Protein Resource (UniProt)(GRCh38).
                                "{3}" stands for Ensembl id from Gencode version 31(GRCh38); 
                                "{4}" stands for Transcript id from Gencode version 31(GRCh38);  
                                "{5}" stands for Ensembl id from Ensembl release 75 version(GRCh37).
                                "{6}" stands for Transcript id from Ensembl release 75 version(GRCh37).
                                "{7}" stands for HGNC id from \x1B[3mHUGO Gene Nomenclature Committee\x1B[23m.
                                "{8}" provides the chromosome or cytogenetic band and extracted from HGNC, Entrez Gene or Ensembl. 
                                "{9}" displays descriptions of a Gene type.
                                "{10}" provides the source of a Gene, like Ens_v31 (extracted from Gencode version31),Retired (Ensembl release-75).

    -o/--output <FILE>          Write output to FILE [stdout] if input is a file. 
                                Default: (date+hour+minute+'_Id_Transform_result.txt'), like "202012241526_Id_Transform_result.txt".
""".format(os.path.basename(sys.argv[0]),str(HEADER_ALL[:2]).rstrip(']'),str(HEADER_ALL[2:10]).strip('[]'),HEADER_ALL[6],HEADER_ALL[7],HEADER_ALL[8],HEADER_ALL[9],HEADER_ALL[10],HEADER_ALL[11],HEADER_ALL[12],HEADER_ALL[13],str(HEADER_ALL[10:]).lstrip('['),HEADER_ALL[3],HEADER_ALL[4],HEADER_ALL[5]))

def version():
    print("""Author: jiabg
Version: 8.0.0
Updated date: 2021-01-15""")


def id_vs_index(index_path):
    '''
    index id with index path
    '''
    global idvsindex,Ensembl_list,Ens_75_index,Ens_v31_index,Gene_info_index,Uniport_info_index
    idvsindex={}
    for root, dirs, files in os.walk(index_path):
        for i in dirs:
            if re.search('Ens_75_from_Ensembl',i):
                idvsindex['Ensembl_75']={}
                idvsindex['Ensembl_75']['index_dir']=os.path.join(root,i)
                ix=open_dir(idvsindex['Ensembl_75']['index_dir'])
                idvsindex['Ensembl_75']['index']=[i for i in ix.schema.names.__self__._subfields.keys()]
            elif re.search('Ens_v31_from_Gencode',i):
                idvsindex['Gencode_v31']={}
                idvsindex['Gencode_v31']['index_dir']=os.path.join(root,i)
                ix=open_dir(idvsindex['Gencode_v31']['index_dir'])
                idvsindex['Gencode_v31']['index']=[i for i in ix.schema.names.__self__._subfields.keys()]
            elif re.search('Gene_info_from_NCBI_and_HGNC',i):
                idvsindex['Gene_info']={}
                idvsindex['Gene_info']['index_dir']=os.path.join(root,i)
                ix=open_dir(idvsindex['Gene_info']['index_dir'])
                idvsindex['Gene_info']['index']=[i for i in ix.schema.names.__self__._subfields.keys()]
            elif re.search('Uniport_info_GRCh38',i):
                idvsindex['Uniport_info']={}
                idvsindex['Uniport_info']['index_dir']=os.path.join(root,i)
                ix=open_dir(idvsindex['Uniport_info']['index_dir'])
                idvsindex['Uniport_info']['index']=[i for i in ix.schema.names.__self__._subfields.keys()]
    Ensembl_list=[idvsindex['Gencode_v31']['index_dir'],idvsindex['Ensembl_75']['index_dir']]
    final_index_header=set(idvsindex['Uniport_info']['index']+idvsindex['Gene_info']['index']+idvsindex['Gencode_v31']['index']+idvsindex['Ensembl_75']['index'])
    for i in ['Joye_id','ENSG_v31','ENSG_75','Symbol_v31','Symbol_75']:
        final_index_header.discard(i)
    return list(final_index_header)

def get_gid_out(joyeid,cindex):
    '''Only for search Joye_id'''
    ix=open_dir(cindex) ##index dir
    searcher=ix.searcher()
    results=searcher.find('Joye_id',joyeid)
    return results

def get_all_out(abc,gstr,cindex):
    '''
    For all search item except Joye_id
    '''
    if abc in ['Ens_v31','Ens_75']:
        for i in Ensembl_list:
            ix=open_dir(i)
            searcher=ix.searcher()
            results=searcher.find(abc,gstr,limit=80)
            abc='Ens_75'
            if len(results)>0:
                break
    elif abc in ['Trans_Ens_v31','Trans_Ens_75']:
        for i in Ensembl_list:
            ix=open_dir(i)
            searcher=ix.searcher()
            results=searcher.find(abc,gstr,limit=80)
            abc='Trans_Ens_75'
            if len(results)>0:
                break
    elif abc in ['ENSG_v31','ENSG_75']:
        for i in Ensembl_list:
            ix=open_dir(i)
            searcher=ix.searcher()
            results=searcher.find(abc,gstr,limit=80)
            abc='ENSG_75'
            if len(results)>0:
                break
    else:
        ix=open_dir(cindex)
        searcher=ix.searcher()
        results=searcher.find(abc,gstr,limit=80)
    return results

def deal_single(input_str):
    '''
    The input is a string 
    '''
    input_str=input_str.strip().strip('"').strip("\'")
    s={}
    if re.search('^\d+$',input_str):
        rs=get_all_out('Gene_id',input_str,idvsindex['Gene_info']['index_dir'])
        s=deal_search_results_all(rs,'Gene_id',input_str,'Accurate')
    elif re.match(r'[OPQ][0-9][A-Z0-9]{3}[0-9]$|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}$', input_str)  and len(input_str) >=6:
        rs=get_all_out('Uniport_id',input_str,idvsindex['Uniport_info']['index_dir'])
        s=deal_search_results_all(rs,'Uniport_id',input_str,'Accurate')
    elif re.match(r'ENSG', input_str):
        if re.search('\.',input_str):
            rs=get_all_out('Ens_v31',input_str,idvsindex['Gencode_v31']['index_dir'])
            s=deal_search_results_all(rs,'Ens_v31',input_str,'Accurate')
        else:
            rs=get_all_out('ENSG_v31',input_str,idvsindex['Gencode_v31']['index_dir'])
            s=deal_search_results_all(rs,'ENSG_v31',input_str,'Accurate')
    elif re.match(r'ENST', input_str):
        if re.search('\.',input_str):
            rs=get_all_out('Trans_Ens_v31',input_str,idvsindex['Gencode_v31']['index_dir'])
        elif re.search('(ENST\d+)',input_str):
            tem=re.search('(ENST\d+)',input_str).group(1)+'??*'
            ix=open_dir(idvsindex['Gencode_v31']['index_dir'])
            searcher=ix.searcher()
            qp=QueryParser('Trans_Ens_v31', schema=ix.schema, termclass=FuzzyTerm)
            q=qp.parse(tem)
            rs=searcher.search(q,limit=None)
        else:
            rs=''
        s=deal_search_results_all(rs,'ENST_v31',input_str,'Accurate')
    elif re.match(r'YP_|XP_|NP_', input_str): #refseq_id  {'YP', 'XP', 'NP'}   transcript  {'XR', 'NR'}
        rs=get_all_out('Refseq_id_pro',input_str,idvsindex['Uniport_info']['index_dir'])
        s=deal_search_results_all(rs,'Refseq_id_pro',input_str,'Accurate')
    elif re.match(r'XR_|NR_|NM_|XM_', input_str):
        rs=get_all_out('Trans_NCBI',input_str,idvsindex['Gene_info']['index_dir'])
        s=deal_search_results_all(rs,'Trans_NCBI',input_str,'Accurate')
    else:
        rs=get_all_out('Symbol',input_str,idvsindex['Gene_info']['index_dir'])
        if len(rs)>=1:
            s=deal_search_results_all(rs,'Symbol',input_str,'Accurate')
        else:
            rs=get_all_out('Alias',input_str,idvsindex['Gene_info']['index_dir'])
            s=deal_search_results_all(rs,'Alias',input_str,'Accurate')
    return s

def deal_file(input_file,ouput_file,out_header,header,print_open):
    '''The input is a file'''
    out=open(ouput_file,'w')
    sep='\t'
    if re.search('csv$',str(input_file)): sep=','
    line_c=1
    filelen=int(os.popen('wc -l {}'.format(input_file)).read().split()[0])
    res_print=set()
    with open(input_file) as f:
        for i in f:
            l=i.strip("\n").split(sep, 1)
            t=l[0].strip()
            if line_c==1 and not re.search('^#',i):
                out.write("%s\t%s\n" % ("Input",sep.join(out_header)))
            if line_c==1 and re.search('^#',i):
                filelen-=1
                if len(l)>1:
                    out.write("%s\t%s\t%s\n" % (t,sep.join(out_header),sep.join(l[1:])))
                else:
                    out.write("%s\t%s\n" % (t,sep.join(out_header)))
            else:
                d=display_results_for_file(deal_single(l[0].strip().strip('"').strip("\'")),l[0],out_header,header,print_open)
                if len(d)>=1:
                    for m in out_header:
                        if not d[m]:
                            d[m]='-'
                        t=t+sep+';'.join(d[m])
                    res_print.add(t)
                    if len(l)>1:
                        out.write("%s\n" % (sep.join([t,l[1]])))
                    else:
                        out.write("%s\n" % (t))
                else:
                    if len(l)>1:
                        out.write("%s\n" % (sep.join([t,sep.join('-'*len(out_header)),l[1]])))
                    else:
                        out.write("%s\n" % (sep.join([t,sep.join('-'*len(out_header))])))
            line_c+=1
    count=0        
    if print_open:
        if len(res_print)>30:
            print('Showing 30 of %d Results for "%s"' % (len(res_print),input_file))
            print("%s\t%s" % ("Input",sep.join(out_header)))
            for j in res_print:
                print(j)
                count+=1
                if count>=30:
                    break
        else:
            print('Showing %d Results for "%s"' % (len(res_print),input_file)) if len(res_print)>1 else print('Showing %d Result for "%s"' % (len(res_print),input_file))
            print("%s\t%s" % ("Input",sep.join(out_header)))
            for j in res_print:
                print(j)
        print('')

def display_results_for_file(result,single,out_header,header,print_open):
    '''
    '''
    s={}
    if len(result)>0:
        for i in result:
            t=single
            for j in out_header:
                t=t+"\t"+result[i][j]
                if j not in s:
                    s[j]=set()
                s[j].add(result[i][j])
    return s

def display_results_for_single(result,single,out_header,header,print_open):
    '''
    Only print result to the screen.
    '''
    fin_header=''
    d=set()
    count=0
    #print(result,out_header)    ##20210118
    if len(result)>0:
        fin_header="Input\t"+'\t'.join(out_header)
        for i in result:
            t=single
            for j in out_header:
                t=t+"\t"+result[i][j]
            d.add(t)
    if len(d)<=0:
        print('Accurate Search: no record for "%s"' % (single))
    elif len(d)>20:
        print('Accurate Search: Showing 20 of %d Results for "%s"' % (len(d),single))
        print(fin_header)
        for i in d:
            print(i)
            count+=1
            if count>20: break
    else:
        print('Accurate Search: Showing %d Results for "%s"' % (len(d),single)) if len(d)>1 else print('Accurate Search: Showing %d Result for "%s"' % (len(d),single))
        print(fin_header)
        for i in d:
            print(i)
    print('')
    return d

def fuzzy_search_for_single(input_str,out_header,header):
    '''
    For fuzzy search...
    '''
    if re.search('^\d+$',input_str):
        input_tem=input_str+'~1'
        d=fuzzy_search_display(input_tem,out_header,header,'Gene_id',idvsindex['Gene_info']['index_dir'],input_str)
    elif re.match(r'[OPQ][0-9][A-Z0-9]{3}[0-9]$|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}$', input_str)  and len(input_str) >=6:
        input_tem=input_str+'~1'
        d=fuzzy_search_display(input_tem,out_header,header,'Uniport_id',idvsindex['Uniport_info']['index_dir'],input_str)
    elif re.match(r'ENSG', input_str): 
        input_tem=input_str+'~1'
        if len(input_str)<14:
            print('---------------------------')
            print('')
            print('Warning: The input "{}" is too short for Fuzzy Searching.'.format(input_str))
            print('')
            sys.exit(0)
        else:
            d=fuzzy_search_display(input_tem,out_header,header,'Ens_v31',idvsindex['Gencode_v31']['index_dir'],input_str)
    elif re.match(r'ENST', input_str):
        if len(input_str)<14:
            print('---------------------------')
            print('')
            print('Warning: The input "{}" is too short for Fuzzy Searching.'.format(input_str))
            print('')
            sys.exit(0)
        else:
            input_tem=input_str+'~1'
            d=fuzzy_search_display(input_tem,out_header,header,'ENST_v31',idvsindex['Gencode_v31']['index_dir'],input_str)
    elif re.match(r'YP_|XP_|NP_', input_str): 
        input_tem=input_str+'~1'
        d=fuzzy_search_display(input_tem,out_header,header,'Refseq_id_pro',idvsindex['Uniport_info']['index_dir'],input_str)
    elif re.match(r'XR_|NR_|NM_|XM_', input_str): 
        input_tem=input_str+'~1'
        d=fuzzy_search_display(input_tem,out_header,header,'Trans_NCBI',idvsindex['Gene_info']['index_dir'],input_str)
    else:
        input_tem=input_str+'~1'
        d=fuzzy_search_display(input_tem,out_header,header,'Symbol',idvsindex['Gene_info']['index_dir'],input_str)
        if len(d)<=0:
            d=fuzzy_search_display(input_tem,out_header,header,'Alias',idvsindex['Gene_info']['index_dir'],input_str)
    print('')
    print('---------------------------')
    print('')
    fres=set()
    count=0
    fhead=''
    if len(d)>0:
        fhead='Input'+"\t"+"\t".join(out_header)
        for i in d:
            t=input_str
            for j in out_header:
               t=t+"\t"+d[i][j]
            fres.add(t)
    if len(fres)<=0:
        print('Fuzzy search: no record for "%s"' % (input_str))
    elif len(fres)>30:
        print('Fuzzy search: Showing 30 of %d Results for "%s"' % (len(fres),input_str))
        print(fhead)
        for m in fres:
            print(m)
            count+=1
            if count==30:
                break
    else:
        print('Fuzzy search: Showing %d Results for "%s"' % (len(fres),input_str)) if len(fres)>1 else print('Fuzzy search: Showing %d Result for "%s"' % (len(fres),input_str))
        print(fhead)
        for m in fres:
            print(m)
    print('')

def fuzzy_search_display(input_tem,out_header,header,index,index_dir,input_str):
    ix = open_dir(index_dir)
    searcher = ix.searcher()
    qp = QueryParser(index, schema=ix.schema)
    qp.add_plugin(qparser.FuzzyTermPlugin())
    q = qp.parse(input_tem)
    res = searcher.search(q,limit=120)
    s=deal_search_results_all(res,index,input_str,'Fuzzy')
    return s

def deal_search_results_all(res,search_index,input_str,fuzzy_or_not):
    d={}
    count=0
    index_l=['Ensembl_75', 'Gencode_v31', 'Gene_info', 'Uniport_info']
    if len(res)<1:
        return d
    for i in res:
        s_count=0
        if search_index=='Symbol' and i['Symbol']==input_str:
            s_count+=1
        elif search_index=='Alias':
            if re.search(';',i['Alias']):
                for m in i['Alias'].strip().split(';'):
                    if m==input_str: s_count+=1
            elif i['Alias']==input_str: s_count+=1
        if search_index=='Symbol' or search_index=='Alias':
            if fuzzy_or_not=='Fuzzy' and s_count:continue
            elif fuzzy_or_not=='Accurate':
                if not s_count: continue
        if i['Joye_id'] not in d:
            d[i['Joye_id']]={}
        if 'Ens_75' in i.keys(): 
            for j in idvsindex['Ensembl_75']['index']:
                if j in d[i['Joye_id']]:
                    d[i['Joye_id']][j]=d[i['Joye_id']][j]+'|'+i[j]
                else:
                    d[i['Joye_id']][j]=i[j]
            if 'Ensembl_75' in index_l: index_l.remove('Ensembl_75')
            for m in index_l:
                ras=get_gid_out(i['Joye_id'],idvsindex[m]['index_dir'])
                if len(ras)>0:
                    for n in ras:
                        for k in idvsindex[m]['index']:
                            if k in d[i['Joye_id']]:
                                d[i['Joye_id']][k]=d[i['Joye_id']][k]+'|'+n[k]
                            else:
                                d[i['Joye_id']][k]=n[k]
                else:
                    for k in idvsindex[m]['index']:
                        d[i['Joye_id']][k]='-'
        elif 'Ens_v31' in i.keys():     ##Ens_v31_index
            for j in idvsindex['Gencode_v31']['index']:
                if j in d[i['Joye_id']]:
                    d[i['Joye_id']][j]=d[i['Joye_id']][j]+'|'+i[j]
                else:
                    d[i['Joye_id']][j]=i[j]
            if 'Gencode_v31' in index_l: index_l.remove('Gencode_v31')
            for m in index_l:
                ras=get_gid_out(i['Joye_id'],idvsindex[m]['index_dir'])
                if len(ras)>0:
                    for n in ras:
                        for k in idvsindex[m]['index']:
                            if k in d[i['Joye_id']]:
                                if d[i['Joye_id']][k]!=n[k]:
                                    d[i['Joye_id']][k]=d[i['Joye_id']][k]+'|'+n[k]
                            else:
                                d[i['Joye_id']][k]=n[k]
                else:
                    for k in idvsindex[m]['index']:
                        d[i['Joye_id']][k]='-'
        elif 'Trans_NCBI' in i.keys():  ##Gene_info_index
            for j in idvsindex['Gene_info']['index']:
                if j in d[i['Joye_id']]:
                    if d[i['Joye_id']][k]!=n[k]:
                        d[i['Joye_id']][j]=d[i['Joye_id']][j]+'|'+i[j]
                else:
                    d[i['Joye_id']][j]=i[j]
            if 'Gene_info' in index_l: index_l.remove('Gene_info')
            for m in index_l:
                ras=get_gid_out(i['Joye_id'],idvsindex[m]['index_dir'])
                if len(ras)>0:
                    for n in ras:
                        for k in idvsindex[m]['index']:
                            if k in d[i['Joye_id']]:
                                if d[i['Joye_id']][k]!=n[k]:
                                    d[i['Joye_id']][k]=d[i['Joye_id']][k]+'|'+n[k]
                            else:
                                d[i['Joye_id']][k]=n[k]
                else:
                    for k in idvsindex[m]['index']:
                        d[i['Joye_id']][k]='-'
        elif 'Uniport_id' in i.keys():  ##Uniport_info_index
            for j in idvsindex['Uniport_info']['index']:
                if j in d[i['Joye_id']]:
                    d[i['Joye_id']][j]=d[i['Joye_id']][j]+'|'+i[j]
                else:
                    d[i['Joye_id']][j]=i[j]
            if 'Uniport_info' in index_l: index_l.remove('Uniport_info')
            for m in index_l:
                ras=get_gid_out(i['Joye_id'],idvsindex[m]['index_dir'])
                if len(ras)>0:
                    for n in ras:
                        for k in idvsindex[m]['index']:
                            if k in d[i['Joye_id']]:
                                if d[i['Joye_id']][k]!=n[k]:
                                    d[i['Joye_id']][k]=d[i['Joye_id']][k]+'|'+n[k]
                            else:
                                d[i['Joye_id']][k]=n[k]
                else:
                    for k in idvsindex[m]['index']:
                        d[i['Joye_id']][k]='-'
        count+=1
        if d[i['Joye_id']]['Symbol_v31']!='-' and d[i['Joye_id']]['Symbol_v31']!=d[i['Joye_id']]['Symbol']:
            tlias=d[i['Joye_id']]['Alias'].split(';')
            tlias.append(d[i['Joye_id']]['Symbol'])
            tlias=set(tlias)
            tlias.discard('')
            tlias.discard('-')
            tlias=';'.join(tlias)
            if not tlias:
                tlias='-'
            d[i['Joye_id']]['Alias']=tlias
            d[i['Joye_id']]['Symbol']=d[i['Joye_id']]['Symbol_v31']
    return d

def main():
    global HEADER,HEADER_ALL,FIN_ENS,header_all_index
    header_all_index=id_vs_index('/home/user/Pipeline/202012_Gene_re/20210118_demo_version5')
    HEADER_ALL=['Gene_id','Symbol','Alias','Trans_NCBI','Refseq_id_pro','Uniport_id','Ens_v31','Trans_Ens_v31','Ens_75','Trans_Ens_75', 'HGNC', 'Location' ,'Genetype','Ens_version'] ##
    try:
        opts, args=getopt(sys.argv[1:], "hvi:o:c:ap", ["help", "version","all","input=","custom=","output=","print"])
        if args != []:
            usage()
            print("Unexpected arguments:", args)
            sys.exit(0)
        if opts == []:
            usage()
            print("No argument found!")
            sys.exit(0)
    except GetoptError:
        usage()
        print("Error: unexpected!")
        sys.exit(0)
    ## Set Default values
    single="150"
    out_file=time.strftime("%Y%m%d%H%M", time.localtime())+'_Id_Transform_result.txt'
    out_header=['Symbol']
    FIN_ENS=''   ## ENSG  ENST final index 
    print_open=False
    for opt, arg in opts:
        if opt in ("-h", "help"):
            usage()
            sys.exit(0)
        elif opt in ("-v", "version"):
            version()
            sys.exit(0)
        elif opt in ["-a", "all"]:
            out_header=header_all_index[:]
        elif opt in ["-c", "--custom"]:
            out_header=arg.strip().split(',')
            temh=[]
            for i in out_header:
                if i not in header_all_index:
                    temh.append(i)
            if len(temh):
                print("Error: {} not in {}".format(temh, header_all_index))
                sys.exit(0)
        elif opt in ["-p", "print"]:
            print_open=True
        elif opt in ["-i", "--input"]:
            single=arg
        elif opt in ["-o", "--output"]:
            out_file=arg
    # main pipe
    if os.access(single, os.R_OK):
        deal_file(single,out_file,out_header,header_all_index,print_open)
    else:
        display_results_for_single(deal_single(single),single,out_header,header_all_index,print_open)
        fuzzy_search_for_single(single,out_header,header_all_index)


########################################
## Main
if __name__ == "__main__":
    main()

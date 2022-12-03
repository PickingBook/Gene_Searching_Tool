# Gene_Searching_Tool
*******
*******

### **v1.0.0**

  A flexible demo tool for searching Gene Id, Gene Symbol, Ensembl Id , Uniport Id , Refseq Id of protein and Transcript of Gencode.v31(GRCh38) and Ensembl.75(GRCh37).

******

#### Usage

```shell
  python3 path/{0} -i <Search id> [-c Symbol,Uniport_id] [-a -o <FILE>]
```

Option:

    -a/--all                        Write|Print all ids to output. 

    -p/--print                      Print some results to the screan when input is a file.

    -v/--version                    Version.

    -h/--help                       Print this help menu.


Parameter:

    -i/--input  <INT|STR|FILE>      Search id or file, the file should be separated by "\\t" or "," . Default: [150].

                                    File will be processed by column. The "xlsx" format is not supported for this version, and will be supported next.

    -c/--custom <STR>               Write|Print the specified id(s) (Symbol,Uniport_id,...). 
                                                                                                                                                                                                                                                                                                                                                                                                                     -o/--output <FILE>             Write output to FILE [stdout] if input is a file. 
                                                                                                                                                                                                                                            Default: (date+hour+minute+'_Id_Transform_result.txt'), like "202012241526_Id_Transform_result.txt".



**********

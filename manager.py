import pandas as pd
import subprocess
import sys
import os


class dataset:
    def __init__(self, path, gene_col, fc_col, title=None, pd_kwargs={}, filters=None, col_funs=[], col_names=[], pmid=None):
        self.df = None
        self.path = path
        self.pd_kwargs = pd_kwargs
        self.gene_col = gene_col
        self.fc_col = fc_col
        self.col_funs = col_funs
        self.col_names = col_names
        self.title = title
        self.pmid = pmid
        self.filters = filters

        self.load_df()
        self.upper_gene_col()
        self.filter_data()
        self.drop_col()
        self.apply_cols()
        self.rename_cols()
        self.col_names=[]
        self.df = self.df.rename(columns={self.fc_col: self.fc_col+"_"+self.title})

    def rename_cols(self):
        if type(self.col_names) is not list:
            self.col_names = [self.col_names]
        self.col_names.append((self.gene_col, "Gene ID"))
        #print(self.col_names)
        for fun in self.col_names:
            self.df = self.df.rename(
                columns={old: new for old, new in self.col_names})
        for col_name in self.col_names:
            for key, value in self.__dict__.items():
                if type(value) == str:
                    if col_name[0] == value:
                        setattr(self, key, col_name[1])

    def apply_cols(self):
        if self.col_funs == []:
            return
        if self.col_funs is not list:
            self.col_funs = [self.col_funs]
        for fun in self.col_funs:
            self.df[fun[1]] = self.df[fun[1]].apply(fun[0])

    def at_index(self, df, index):
        if index is int:
            return df.iloc[index]
        else:
            return df[index]

    def load_df(self):
        if self.path.endswith(".xlsx"):
            self.df = pd.read_excel(self.path, **self.pd_kwargs,engine="openpyxl")
        else:
            self.df = pd.read_csv(self.path, **self.pd_kwargs)

    def upper_gene_col(self):
        self.df[self.gene_col] = self.df[self.gene_col].apply(lambda x: str(x).upper())

    def drop_col(self):
        # TODO: specify which columns to keep
        self.df = self.df[[self.gene_col, self.fc_col]]

    def auto_filter(self):
        if "FDR" in self.df.columns:
            return (lambda x : float(x) < 0.05,"FDR")
        if "pval" in self.df.columns:
            return (lambda x : float(x) < 0.05,"pval")
        return None

    def filter_data(self):
        if self.filters == None:
            self.filters = self.auto_filter()
        if self.filters is not list:
            self.filters = [self.filters]
        self.filters = list(filter(None,self.filters))
        for filtr in self.filters:
            def filter_row(row):
                return filtr[0](self.at_index(row, filtr[1]))
            try:
                self.df = self.df[self.df.apply(filter_row, axis=1)]
            except Exception as e:
                print(e)
                print(self.title)
                pass


def merge_datasets(datasets):
    merged = pd.merge(
        left=datasets[0].df, right=datasets[1].df,
        left_on=["Gene ID", "Original Gene ID"], right_on=["Gene ID","Original Gene ID"],how="outer")
            #left_on="Gene ID", right_on="Gene ID",how="outer")
    for i in range(2, len(datasets), 1):
        merged = pd.merge(
            left=merged, right= datasets[i].df,
            left_on=["Gene ID", "Original Gene ID"], right_on=["Gene ID","Original Gene ID"],how="outer")
            #left_on="Gene ID", right_on="Gene ID",how="outer")
    return merged

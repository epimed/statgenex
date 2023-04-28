from src.statgenex import Analysis
from src.statgenex.service import FormatService, FigureService, FileService, DataLoader, IndexReducer
from src.statgenex.stats import BenjaminiHochberg
import pandas as pd
import numpy as np
from scipy.stats import f_oneway, kruskal
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import warnings
import openpyxl

# ==============================

class Anova(Analysis):
    
    def __init__(self,
                 project,
                 dataset_name,
                 group_names,
                 features,
                 **kwargs
                 ):
        super().__init__()
        self.project = project
        self.dataset_name = dataset_name
        self.group_names = group_names
        self.features = features
        
        self.generate_plots = True
        self.generate_pvalues = True
        self.figsize = None
        self.figwidth_scale = 1.2
        self.regular = 16
        self.show_title = True
        self.show_pval_anova = True
        self.show_pval_kw = True
        self.boxplot_options = FigureService.create_boxplot_options()
        
        self.significance = {'pval_anova': 0.05, 'pval_kw': 0.05, 'fdr_anova': 0.05, 'fdr_kw': 0.05}
        
        for k, v in kwargs.items():
            setattr(self, k, v)
    
        self.results = pd.DataFrame(index=self.features)
        self.results.index.name = 'gene'
        self.description = pd.DataFrame()
        self.description.index.name = 'group_name' 
        self.sample_sizes = pd.DataFrame(index=self.features)
        self.aov_data_dict = dict()
    
    @property 
    def name(self):
        return "Anova"
    
    @property
    def results_dir(self):
        self.local_dir = f"{FormatService.today()}_{self.name}/"
        return self.project.results_dir + self.local_dir
    
    def perform(self):
        dataset = self.project.datasets[self.dataset_name]
        expression_data = self._generate_expression_data()
        self.available_group_names = [gn for gn in self.group_names if gn in dataset.groups.keys()]
        for group_name in self.available_group_names:
            self.description.loc[group_name, 'dataset_name'] = self.dataset_name
            self.description.loc[group_name, 'sample_size'] = len(dataset.groups[group_name].samples)
        self._calculate_anova(dataset, expression_data)
        self._calculate_fdr()
        self._calculate_significance()
        if self.generate_plots:
            FileService.create_folder(self.results_dir)
            self._generate_boxplots()
        if self.generate_pvalues:
            self.save_results()
    
    def _generate_boxplots(self):
        n_groups = len(self.available_group_names)
        figwidth = self.figwidth_scale*n_groups
        figsize = (figwidth, 4) if self.figsize is None else self.figsize
        self.pdf_filename = self.results_dir + f"Anova_boxplots_{self.dataset_name}_{len(self.features)}_genes_{len(self.available_group_names)}_groups.pdf"
        with PdfPages(self.pdf_filename) as pdf:
            for feature in self.features:
                if feature in self.aov_data_dict.keys():
                    aov_data = self.aov_data_dict[feature]
                    fig, ax = plt.subplots(figsize=figsize)
                    ax.boxplot(aov_data, **self.boxplot_options)                    
                    self._add_annotations(ax, feature)
                    pdf.savefig(fig, bbox_inches='tight', orientation='landscape')
                    plt.close('all')
    
    def _add_annotations(self, ax, feature):
        font = FigureService.create_arial_narrow_font()
        regular, medium, small, tiny = FigureService.create_font_sizes(regular=self.regular)
        
        title = feature + ' - ' + self.dataset_name
        if self.show_pval_anova:
            pval = self.results.loc[feature, 'pval_anova']
            pAnovaText = 'ANOVA p-value = ' + '{:.1e}'.format(pval) + ' ' + FigureService.get_significance_symbol(pval)
            pAnovaText = pAnovaText.strip()
            title = title + '\n' + pAnovaText
        
        if self.show_pval_kw:
            pval = self.results.loc[feature, 'pval_kw']
            pAnovaText = 'Kruskal-Wallis p-value = ' + '{:.1e}'.format(pval) + ' ' + FigureService.get_significance_symbol(pval)
            pAnovaText = pAnovaText.strip()
            title = title + '\n' + pAnovaText
        
        if self.show_title:
            ax.set_title(title, fontsize=regular, **font)
            
        xticks = [i+1 for i in range(len(self.available_group_names))]
        ax.set_xticks(xticks)
        
        xticklabels = []
        for group_name in self.available_group_names:
            # n_samples = len(self.project.datasets[self.dataset_name].groups[group_name].samples)
            n_samples = self.sample_sizes.loc[feature, group_name]
            xticklabel = group_name + '\n' + '(n=' + '{:.0f}'.format(n_samples) + ')' 
            xticklabels.append(xticklabel)
        ax.set_xticklabels(xticklabels, fontsize=medium, **font)
        
        ax.set_ylabel('Expression', fontsize=regular, **font)
        ax.tick_params(axis='y', labelsize=tiny)    
    
    def _calculate_significance(self):
        query = True
        for k, v in self.significance.items():
            query = query & (self.results[k]<v)
        self.results.loc[query, 'significant'] = 1
        self.results.loc[~query, 'significant'] = 0
        
    def _calculate_fdr(self):
        bh = BenjaminiHochberg()
        for test_name in ['anova', 'kw']:
            p_vals = self.results['pval_' + test_name]
            fdr = bh.perform(p_vals)
            self.results['fdr_' + test_name] = fdr
        
    def _calculate_anova(self, dataset, expression_data):
        for feature in expression_data.index:
            aov_data = []
            for group_name in self.available_group_names:
                group_samples = dataset.groups[group_name].samples
                not_na_values = expression_data.loc[feature, group_samples].dropna()
                self.sample_sizes.loc[feature, group_name] = len(not_na_values)
                aov_data.append(list(not_na_values))
            self.aov_data_dict[feature] = aov_data
            try:
                with warnings.catch_warnings(record=True):
                    warnings.simplefilter("always")
                    f_aov, pval_aov = f_oneway(*aov_data)
                    f_kw, pval_kw = kruskal(*aov_data)
                    self.results.loc[feature, 'pval_anova'] = pval_aov
                    self.results.loc[feature, 'pval_kw'] = pval_kw
            except:
                pass
        
    def _generate_expression_data(self):
        dataset = self.project.datasets[self.dataset_name]
        data_loader = DataLoader(dataset.data_dir + dataset.data_filename)
        data_loader.load()
        reducer = IndexReducer(data=data_loader.data, features=self.features)
        return reducer.transform()

    def save_results(self):
        FileService.create_folder(self.results_dir)
        # self.results.to_csv(self.results_dir + 'Anova_pvalues.csv', sep=';', index=True)
        # self.description.to_csv(self.results_dir + 'Anova_sample_sizes.csv', sep=';', index=True)
        output_prefix = f"Anova_results_{self.dataset_name}_{len(self.features)}_genes_{len(self.available_group_names)}_groups"
        with pd.ExcelWriter(self.results_dir + output_prefix + '.xlsx', engine='openpyxl') as writer:
            self.results.to_excel(writer, sheet_name='p-values')
            self.description.to_excel(writer, sheet_name='sample_sizes')
            significance = pd.DataFrame()
            significance.index.name = 'pval_type'
            for k, v in self.significance.items():
                significance.loc[k, 'threshold'] = v
            significance.to_excel(writer, sheet_name='significance')
        
    def __repr__(self):
        return (f"{self.__class__.__name__} ["
                f"name = {self.name}, "
                f"project_name = {self.project.name}, "
                f"dataset_name = {self.dataset_name}, "
                f"group_names = {self.group_names}, "
                f"features = {self.features}"
                f"]")

# ==============================

from src.statgenex.service import FormatService, FileService, DataLoader
import json
from abc import ABC, abstractmethod

# ==============================

class Persistent(ABC):
    """
    Persistent interface:
    Allow to persist an object as a JSON file.
    """

    @abstractmethod
    def dump(self) -> None:
        """Dump the object as a JSON file"""
        ...
    
    @abstractmethod    
    def restore(self) -> None:
        """Restore the object from a JSON file"""
        ...

# ==============================

class Entity:
    """Abstract Entity"""
    
    def as_dict(self):
        as_dict = dict()
        for k, v in self.__dict__.items():
            if isinstance(v, dict):
                sub_dict = dict()
                for ki, vi in v.items():
                    sub_dict[ki] = vi.as_dict()
                as_dict[k] = sub_dict
            else:
                as_dict[k] = v
        return as_dict
    
    def __repr__(self):
        attributes = ", ".join("=".join((str(k), str(v))) for k, v in self.__dict__.items())
        return (f"{self.__class__.__name__} [{attributes}]")

# ==============================

class Project(Entity, Persistent):
    def __init__(self, name, root_dir, **kwargs):
        self.name = FormatService.normalize(name)
        self.root_dir = root_dir
        self.project_dir = self.root_dir + name + '/'
        self.data_dir = self.project_dir + 'data/'
        self.results_dir = self.project_dir + 'results/'
        self.json_file = self.project_dir + 'project.json'
        self.datasets = dict()
        for k, v in kwargs.items():
            if k.endswith('_dir'):
                v = FormatService.normalize_directory_path(v)
            setattr(self, k, v)
    
    def add_dataset(self, dataset):
        dataset.data_dir = self.data_dir
        self.datasets[dataset.name] = dataset
          
    def print_summary(self):
        print('Project', self.name, self.project_dir)
        for dk, dv in self.datasets.items():
            if len(dv.groups)>0:
                group_repr = ", ".join((f"{k} ({len(v.samples)})" for k, v in dv.groups.items()))
                print(f"Dataset {dk}: {group_repr}")
            else:
                print(f"Dataset {dk}: no groups defined")   
   
    def dump(self) -> None:
        folders = []
        for k, v in self.__dict__.items():
            if k.endswith('_dir'):
                folders.append(v)
        for folder in folders:
            FileService.create_folder(folder)
        with open(self.json_file, 'w', encoding='utf-8') as f:
            json.dump(self.as_dict(), f, ensure_ascii=True, indent=4)  
     
    def restore(self):
        """Load project from file project.json"""
        with open(self.json_file, 'r', encoding='utf-8') as f:
            project_dict = json.load(f)   
        for k, v in project_dict.items():
            if (k!='datasets'):
                setattr(self, k, v)
        self.datasets.clear()
        for dataset_dict in project_dict['datasets'].values():
            dataset = Dataset(**dataset_dict)
            dataset.groups = dict()
            for group_dict in dataset_dict['groups'].values():
                group = Group(**group_dict)
                dataset.add_group(group)
            self.add_dataset(dataset)       
    
    def __repr__(self):
        return (f"{self.__class__.__name__} ["
                f"name = {self.name}, "
                f"root_dir = {self.root_dir}, "
                f"project_dir = {self.project_dir}, "
                f"data_dir = {self.data_dir}, "
                f"results_dir = {self.results_dir}"
                f"]")

# ==============================

class Dataset(Entity):
    def __init__(self, name: str, **kwargs):
        self.name = name
        self.data_dir = None
        self.data_filename = None
        self.data_ext = 'csv'
        self.data_sep = ';'
        self.expgroup_filename = None
        self.expgroup_ext = 'csv'
        self.expgroup_sep = ';'
        self.groups = dict()
        for k, v in kwargs.items():
            setattr(self, k, v)
        
    def add_group(self, group: 'Group') -> None:
        self.groups[group.name] = group
    
    def add_groups(self, groups: dict[str, 'Group']) -> None:
        self.groups = groups
    
    def generate_groups(self, categorical_filters=None, quantitative_filters=None, expression_filters=None, secondary_filters=None):
        expgroup_loader = DataLoader(filename=self.data_dir + self.expgroup_filename, ext=self.expgroup_ext, sep=self.expgroup_sep)
        expgroup_loader.load()
        expgroup = expgroup_loader.data
        data_loader = DataLoader(filename=self.data_dir + self.data_filename, ext=self.expgroup_ext, sep=self.expgroup_sep)
        data_loader.load()
        expression_data = data_loader.data
        expression_data = expression_data.dropna(axis=1, how='all')
        expression_data = expression_data.dropna(axis=0, how='all')
        expression_data = expression_data.drop_duplicates()
        common_samples = set(expgroup.index).intersection(set(expression_data.columns))
        expression_data = expression_data[common_samples]
        expgroup = expgroup.loc[common_samples]
        if categorical_filters is not None:
            self._generate_categorical_groups(expgroup, categorical_filters)
        if quantitative_filters is not None:
            self._generate_quantitative_groups(expgroup, quantitative_filters)
        if expression_filters is not None:
            self._generate_expression_groups(expression_data, expression_filters)
        if secondary_filters is not None:
            self._generate_secondary_groups(secondary_filters)
    
    def _generate_expression_groups(self, expression_data, expression_filters):
        for group_name, expression_filter in expression_filters.items():
            ref_group_name = expression_filter['ref_group']
            ref_group_samples = self.groups[ref_group_name].samples
            common_samples = set(expression_data.columns).intersection(set(ref_group_samples))
            gene_name = expression_filter['gene']
            if gene_name in expression_data.index:
                if (expression_filter['threshold_type']=='median'):
                    threshold = expression_data.loc[gene_name, common_samples].median()
                    query = None
                    if expression_filter['class']=='low':
                        query = (expression_data.loc[gene_name, common_samples]<=threshold)
                    if expression_filter['class']=='high':
                        query = (expression_data.loc[gene_name, common_samples]>threshold)
                    group_samples = list(query.loc[query].index)
                    self.add_group(Group(name=group_name, samples=group_samples))
        
    def _generate_quantitative_groups(self, expgroup, quantitative_filters):
        for group_name, list_filters in quantitative_filters.items():
            quantitative_query = self._get_quantitative_query(expgroup, list_filters)
            group_samples = list(expgroup.loc[quantitative_query].index)
            self.add_group(Group(name=group_name, samples=group_samples))
    
    def _get_quantitative_query(self, expgroup, list_filters):
        query_and = True
        for filter_element in list_filters:
            for colname, colvalues in filter_element.items():
                condition_min = (expgroup[colname]>=min(colvalues))
                condition_max = (expgroup[colname]<max(colvalues))
                query_and = query_and & ( condition_min & condition_max )
        return query_and
            
    def _generate_categorical_groups(self, expgroup, categorical_filters):
        for group_name, list_filters in categorical_filters.items():
            categorical_query = self._get_categorical_query(expgroup, list_filters)
            group_samples = list(expgroup.loc[categorical_query].index)
            self.add_group(Group(name=group_name, samples=group_samples))
    
    def _get_categorical_query(self, expgroup, list_filters):
        query_and = True
        for filter_element in list_filters:
            for colname, colvalues in filter_element.items():
                if not isinstance(colvalues, list):
                    colvalues = [colvalues]
                query_and = query_and & (expgroup[colname].isin(colvalues))
        return query_and

    def _generate_secondary_groups(self, secondary_filters):
        for secondary_group_name, list_primary_groups in secondary_filters.items():
            common_samples = set()
            for i, primary_group_name in enumerate(list_primary_groups):
                if (i==0):
                    common_samples = common_samples.union(set(self.groups[primary_group_name].samples))
                else:
                    common_samples = common_samples.intersection(set(self.groups[primary_group_name].samples))
            self.add_group(Group(name=secondary_group_name, samples=sorted(list(common_samples))))    
    
    
    
    def __repr__(self):
        group_names = list(self.groups.keys())
        return (f"{self.__class__.__name__} ["
                f"name = {self.name}, "
                f"data_filename = {self.data_filename}, "
                f"expgroup_filename = {self.expgroup_filename}, "
                f"groups = {group_names}"
                f"]")
    
# ==============================

class Group(Entity):
    """Group of samples"""
    
    def __init__(self, name, **kwargs):
        self.name = FormatService.normalize(name)
        self.samples = []
        for k, v in kwargs.items():
            setattr(self, k, v)
    
    def __repr__(self):
        samples_repr = str(self.samples)
        if len(self.samples)>3:
            samples_repr = str(self.samples[0:3]) + '...'
        return (f"{self.__class__.__name__} [name={self.name}, n={len(self.samples)}, samples={samples_repr}]")
        
# ==============================

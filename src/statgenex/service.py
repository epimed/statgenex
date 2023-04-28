from abc import ABC, abstractmethod
import pandas as pd
from datetime import date, datetime
import matplotlib as mpl
import matplotlib.cm as cm
import matplotlib.colors as clr
import warnings
import openpyxl
import os

# ============================== 

class FormatService:
    
    @classmethod
    def normalize(cls, text):
        return '_'.join(text.strip().split())
    
    @classmethod
    def normalize_lower(cls, text):
        return '_'.join(text.strip().lower().split())
    
    @classmethod
    def today(cls): 
        return date.today().strftime("%Y.%m.%d")
    
    @classmethod
    def now(cls): 
        return datetime.now().strftime("%Y.%m.%d-%H.%M.%S")
    
    @classmethod
    def normalize_directory_path(cls, text):
        if text[-1]!='/':
            text = text + '/'
        return text
        
    

# ==============================    
    
class FileService:
    
    @classmethod
    def create_folder(cls, folder):
        if not os.path.exists(folder):
            os.makedirs(folder)
            
    
# ==============================    
    
class FigureService:
    
    @classmethod
    def get_significance_symbol(cls, pvalue, oneStar=0.05, twoStars=0.01, threeStars=0.001):
        symbol = ''
        if (pvalue<=oneStar):
            symbol = '*'
        if (pvalue<=twoStars):
            symbol = '**'
        if (pvalue<=threeStars):
            symbol = '***'
        return symbol

    @classmethod
    def create_font_sizes(cls, regular=20):
        medium = 0.8 * regular
        small = 0.7 * regular
        tiny = 0.6 * regular
        return regular, medium, small, tiny
    
    @classmethod
    def create_arial_narrow_font(cls):
        mpl.rcParams['font.family'] = 'Arial'
        mpl.rc('font',family='Arial')
        return {'fontname':'Arial', 'stretch' : 'condensed'}
    
    @classmethod
    def create_arial_font(cls):
        mpl.rcParams['font.family'] = 'Arial'
        mpl.rc('font',family='Arial')
        return {'fontname':'Arial'}
    
    @classmethod
    def save_fig_with_resolution(cls, fig, output_dir, file_prefix, dpi=100, ext='png'):
        FileService.create_folder(output_dir)
        filename = output_dir + file_prefix + '.' + ext
        fig.savefig(filename, dpi=dpi, format=ext, bbox_inches='tight', orientation='portrait')
    
    @classmethod    
    def extract_colors_from_colormap(cls, n=10, colormap='jet'):
        cmap = cm.get_cmap(colormap)
        norm = mpl.colors.Normalize(vmin=0, vmax=n-1) 
        return [cmap(norm(ind)) for ind in range(n)] 
    
    @classmethod
    def generate_colors_from_colormap(cls, values, colormap='jet', vmin=None, vmax=None):
        cmap = cm.get_cmap(colormap)
        if (vmin is None):
            vmin = min(values)
        if (vmax is None):
            vmax = max(values)
        norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax) 
        return [cmap(norm(v)) for v in values]  
    
    @classmethod
    def create_custom_colormap(cls, palette='white', n_segments=256, list_colors=None):
        if list_colors is None:
            if palette=='black':
                list_colors = ['cyan', 'royalblue', 'black', 'crimson', 'pink']
            else:
                list_colors = ['royalblue', 'cyan', 'azure', 'whitesmoke', 'lavenderblush', 'pink', 'crimson']
        return clr.LinearSegmentedColormap.from_list('custom', list_colors, N=n_segments)
    
    @classmethod
    def create_boxplot_options(cls):
        '''
        boxprops = dict(linestyle='-', linewidth=0.75, color='black')
        flierprops = dict(marker='o', markersize=4, markeredgewidth=0.3, markeredgecolor='black')
        medianprops = dict(linestyle='-', linewidth=0.75, color='black')
        meanprops = dict(linestyle='-', linewidth=1.5, color='black')
        capprops = dict(color='black', linewidth=0.75)
        whiskerprops = dict(linestyle='--', linewidth=0.75, color='black')            
        boxplot_options = {
            'flierprops': flierprops, 
            'medianprops': medianprops, 
            'meanprops': meanprops, 
            'meanline': True, 
            'showmeans': True,
            'boxprops': boxprops, 
            'capprops': capprops, 
            'whiskerprops': whiskerprops,
            'patch_artist': True,
            'widths': 0.5
            }
        '''    
        boxprops = dict(linestyle='-', linewidth=0.75, color='black', facecolor='silver', alpha=0.5)
        flierprops = dict(marker='o', markersize=2, markeredgewidth=0.3, markeredgecolor='black')
        medianprops = dict(linestyle='-', linewidth=0.75, color='black')
        meanprops = dict(linestyle='-', linewidth=1.5, color='black')
        capprops = dict(color='black', linewidth=0.75)
        whiskerprops = dict(linestyle='--', linewidth=0.75, color='black')            
        boxplot_options = {
            'flierprops': flierprops, 
            'showfliers': False,
            'medianprops': medianprops, 
            'meanprops': meanprops, 
            'meanline': False, 
            'showmeans': False,
            'boxprops': boxprops, 
            'capprops': capprops, 
            'whiskerprops': whiskerprops,
            'patch_artist': True,
            'widths': 0.6,
            'whis': 100,
            'autorange': True,
            }
        return boxplot_options
    

# ==============================

class Loader(ABC):
    """Data loader interface"""
    
    def __init__(self):
        ...
           
    def load(self):
        ...
        
# ==============================

class DataLoader(Loader):
    """Load data from a file into a standard pandas DataFrame"""
    
    def __init__(self, filename, ext='csv', sep=';', sheet_name=0):
        super().__init__()
        self.filename = filename
        self.ext = ext
        self.sep = sep  
        self.sheet_name = sheet_name  
        self.data = None  

    def load(self):
        if self.ext=='excel':
            with warnings.catch_warnings(record=True):
                warnings.simplefilter("always")
                self.data = pd.read_excel(self.filename, engine="openpyxl", sheet_name=self.sheet_name, index_col=0)
        else:
            self.data = pd.read_csv(self.filename, sep=self.sep, index_col=0)

# ==============================

class Transformer(ABC):
    """Interface Data Transformer"""
 
    @abstractmethod
    def transform(self) -> pd.DataFrame:
        ...

# ==============================

class ColumnReducer(Transformer):
    """Reduce columns to a given list"""

    def __init__(self, data: pd.DataFrame, features: list):
        self.data = data
        self.features = features
    
    def transform(self) -> pd.DataFrame:
        if self.features:
            available_features = sorted(list(set(self.data.columns).intersection(set(self.features))))
            return self.data[available_features]
        return self.data
    
# ==============================

class IndexReducer(Transformer):
    """Reduce index to a given list"""
    
    def __init__(self, data: pd.DataFrame, features: list):
        self.data = data
        self.features = features
    
    def transform(self) -> pd.DataFrame:
        if self.features:
            available_features = sorted(list(set(self.data.index).intersection(set(self.features))))
            return self.data.loc[available_features]
        return self.data

# ==============================

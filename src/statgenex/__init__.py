from abc import ABC, abstractmethod

# ==============================       

class Analysis(ABC):
    
    @property
    @abstractmethod
    def name(self) -> str:
        ...
    
    @abstractmethod
    def perform(self):
        ...
    
    @property
    @abstractmethod
    def results_dir(self):
        ...
        
    @abstractmethod
    def save_results(self) -> dict:
        ...

# ============================== 

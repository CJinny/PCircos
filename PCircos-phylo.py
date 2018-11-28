import sys
import plotly.graph_objs as go
from plotly.offline import download_plotlyjs, init_notebook_mode, plot, iplot
init_notebook_mode(connected=True)
import numpy as np
import pandas as pd
import maths
from Complex import Complex
from Bio import Phylo
from time import time
from collections import OrderedDict
from fig import Figure

class circularTree(Figure):
    def __init__(self, *args, **kwargs):
        ## TEMP
        assert 'path' in kwargs
        self.path = kwargs.pop('path')
        if 'format' in kwargs:
            self.format = kwargs.pop('format')
        else:
            self.format = 'phyloxml' 
        
        if 'cladogram' in kwargs:
            self.cladogram = kwargs.pop('cladogram')
            assert isinstance(self.cladogram, bool)
        else:
            self.cladogram = False

        if 'degreerange' in kwargs:
            self.degreerange = kwargs.pop('degreerange')
        else:
            self.degreerange = [0, 360]

        if 'start_leaf' in kwargs:
            self.start_leaf = kwargs.pop('start_leaf')
            assert self.start_leaf in ['first', 'last']
        else:
            self.start_leaf = 'first'
        
        if 'npoints' in kwargs:
            self.npoits = kwargs.pop('npoints')
        else:
            self.npoints = 100
        
        if 'decimals' in kwargs:
            self.decimals = kwargs.pop('decimals')
        else:
            self.decimals = 6

        self.tree = Phylo.read(self.path, self.format)
        self.nleaves = self.tree.count_terminals()
        self.nnodes = len(self.tree.depths())
        self.nroots = self.nnodes - self.nleaves
        self.leaves = self.tree.get_terminals()
        self.roots = self.tree.get_nonterminals()

    def decimals_control(self, dictionary=None, k=None, v=None):
        if dictionary is not None:
            k = list(dictionary.keys())
            v = list(dictionary.values())
        v = np.around(v, decimals=self.decimals)
        
        return OrderedDict(sorted(dict(zip(k,v)).items(), key=lambda t:str(t[0])))

    def check_unit_branch_lengths(self):
        if self.cladogram is True:
            return True
        elif not np.count_nonzero(self.tree.depths().values()):
            return True
        else:
            return False

    def get_node_radius_dict(self):
        
        node_radius = self.tree.depths(unit_branch_lengths=self.check_unit_branch_lengths())
        return self.decimals_control(node_radius)
        
    
    def get_branch_length_dict(self):
        k = list(self.get_node_radius_dict().keys())
        if self.cladogram is True:
            v = np.repeat(1, len(self.get_node_radius_dict()))
        else:
            v =  [*map(lambda x: (0 if x.branch_length is None else x.branch_length), list(self.get_node_radius_dict().keys()))]           
        
        return self.decimals_control(k=k,v=v)

    def get_root_radius_dict(self):
        k = list(self.get_node_radius_dict().keys())
        v1 = np.array(list(self.get_node_radius_dict().values()))
        v0 = np.array(list(self.get_branch_length_dict().values()))
        v = v1 - v0
        return self.decimals_control(k=k,v=v)


    def get_y_dict(self):
        # return dict
        if self.start_leaf == 'last':
            node_y = dict((leaf, k) for k, leaf in enumerate(self.leaves))
        else:
            node_y = dict((leaf, k) for k, leaf in enumerate(reversed(self.leaves)))

        def assign_y(clade):
            for subclade in clade:
                if subclade not in node_y:
                    assign_y(subclade)
            node_y[clade] = 0.5 * (node_y[clade.clades[0]] + node_y[clade.clades[-1]])
        
        if self.tree.root.clades:
            assign_y(self.tree.root)

        return OrderedDict(sorted(node_y.items(), key=lambda t:str(t[0])))

    def y2theta(self, y):
        yvals = self.get_y_dict().values()
        miny, maxy = min(yvals), max(yvals)
        return (np.pi/180) * (self.degreerange[0] + (self.degreerange[1] - self.degreerange[0]) * (y - miny) / float(maxy - miny))
    
    def get_theta_dict(self):
        
        k = list(self.get_y_dict().keys())
        v = [*map(lambda y:self.y2theta(y), list(self.get_y_dict().values()))]
        return OrderedDict(sorted(dict(zip(k, v)).items(), key=lambda t:str(t[0])))


    def get_node_complex_dict(self):
        
        k = list(self.get_node_radius_dict().keys())
        radius = np.array(list(self.get_node_radius_dict().values()))
        theta = np.array(list(self.get_theta_dict().values()))
        Complex = np.zeros(radius.shape, dtype='complex')
        Complex.real = np.sin(theta) * radius
        Complex.imag = np.cos(theta) * radius
        return OrderedDict(sorted(dict(zip(k, Complex)).items(), key=lambda t:str(t[0])))

    def get_root_complex_dict(self):
        k = list(self.get_root_radius_dict().keys())
        radius = np.array(list(self.get_root_radius_dict().values()))
        theta = np.array(list(self.get_theta_dict().values()))
        Complex = np.zeros(radius.shape, dtype='complex')
        Complex.real = np.sin(theta) * radius
        Complex.imag = np.cos(theta) * radius
        return OrderedDict(sorted(dict(zip(k, Complex)).items(), key=lambda t:str(t[0])))

    def get_stem_path(self):
        # from node to its original root
        v0 = np.array(list(self.get_root_complex_dict().values()))
        v1 = np.array(list(self.get_node_complex_dict().values()))
        path = " ".join(np.column_stack((np.repeat('M', len(v0)), v0.real, v0.imag, np.repeat('L',len(v1)), v1.real, v1.imag)).ravel())
        return path

    def get_angular_path(self):
        # resort dict by root radius
        # ONGOING
        # order by root radius
        resort_root_radius = OrderedDict(sorted(self.get_root_radius_dict().items(), key=lambda t: t[1]))
        root_radius_values = np.array(list(resort_root_radius.values()))
        _, ind = np.unique(root_radius_values, return_index=True)

        # sort root theta by root radius values
        resort_theta = OrderedDict(sorted(self.get_theta_dict().items(), key=lambda t: resort_root_radius[t[0]]))
        theta_values = np.array(list(resort_theta.values()))

        radius = np.split(root_radius_values, ind)
        radius.pop(0)
       
        theta_split = np.split(theta_values, ind)
        theta_split.pop(0)

        resolution = 0.5 * self.npoints/np.pi

        radius = [*map(lambda t:t[0], radius)]

        theta_split = [*map(lambda t:np.sort(t), theta_split)]
        theta_range = [*map(lambda t:np.linspace(t[0], t[-1], round(2+(t[-1]-t[0])*resolution)), theta_split)]
        def tocomplex(theta, radius):
            Complex = np.zeros(theta.shape, dtype='complex')
            Complex.real = np.sin(theta) * radius
            Complex.imag = np.cos(theta) * radius
            return Complex
        complex_range = [*map(lambda theta, radius: tocomplex(theta, radius), theta_range, radius)]
        
        # the tree starts from the center of circle, and its complex coord is 0, which means we can remove them
        complex_range.pop(0)
        
        path_list = [*map(lambda x: " ".join(np.column_stack((np.concatenate((np.full(1, 'M'), np.full(len(x)-1, 'L'))), x.real, x.imag)).ravel().tolist()), complex_range)]
        path = " ".join(path_list)

        return path

    def phylo_trace(self):
        node_complex = list(self.get_node_complex_dict().values())
        trace = [go.Scatter(x=[*map(lambda x:x.real, node_complex)], 
                           y=[*map(lambda x:x.imag, node_complex)],
                           mode='markers',
                           hoverinfo='text',
                           text=[*map(lambda t:'{}: {}'.format(t.name, t.branch_length), list(self.get_node_complex_dict().keys()))],
                           marker=dict(color='red', size=10))]
        

        '''
        for i in range(len(node_complex)):
            trace.append(go.Scatter(x=node_complex[i].real,
                                    y=node_complex[i].imag,
                                    marker=dict(color='red', size=10),
                                    opacity=1))
        '''
        
        return trace

    def phylo_layout(self):
        layout = go.Layout(autosize=False,
                           showlegend=False,
                           xaxis=dict(visible=False),
                           yaxis=dict(visible=False),
                           hovermode='closest',
                           width=1500,
                           height=1500)
        layout['shapes'].append(dict(path=self.get_stem_path()))
        layout['shapes'].append(dict(path=self.get_angular_path()))
        return layout

    def phylo_fig(self):
        return go.Figure(data=self.phylo_trace(), layout=self.phylo_layout())

phylo = circularTree(path=sys.argv[1])

t=time()
print(phylo.get_branch_length_dict())
print('\t')
print(phylo.get_node_radius_dict())
print('\t')
print(phylo.get_root_radius_dict())
print('\t')
print(phylo.get_theta_dict())
print('\t')
print(phylo.get_node_complex_dict())
print('\t')
print(phylo.get_root_complex_dict())
print('\t')
print(phylo.get_stem_path())
print('\t')

path = phylo.get_angular_path()
'''
print(radius)
print('\t')
print([*map(lambda t:len(t), theta)])
print('\t')
print(theta)
print('\t')
print(complex_range)
print('\t')
'''

print(path)
print('\t')
print(phylo.phylo_trace())



print('total run time:')
print (time()-t)




fig = phylo.phylo_fig()
plot(fig)
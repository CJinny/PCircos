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
import colorlover as cl

# Dependency:
## Biopython

# important function in biopython Phylo:
# tree.depths()


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
            try:
                assert isinstance(self.cladogram, bool)
            except AssertionError:
                self.cladogram = False
        else:
            self.cladogram = False

        if 'degreerange' in kwargs:
            self.degreerange = kwargs.pop('degreerange')
            try:
                assert isinstance(self.degreerange, list)
                assert self.degreerange[1] >= self.degreerange[0]
            except AssertionError:
                self.degreerange = [90, 90]
        else:
            self.degreerange = [90, 90]

        if 'start_leaf' in kwargs:
            self.start_leaf = kwargs.pop('start_leaf')
            try:
                assert self.start_leaf in ['first', 'last']
            except AssertionError:
                self.start_leaf = 'first'
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

        self.colorscale = [[0.0, 'rgb(10,10,150)'],
                            [0.0009, 'rgb(10,10,150)'],
                            [0.001, 'rgb(214, 47, 38)'],  
                            [0.1, 'rgb(214, 47, 38)'],
                            [0.2, 'rgb(244, 109, 67)'],
                            [0.3, 'rgb(252, 172, 96)'],
                            [0.4, 'rgb(254, 224, 139)'],
                            [0.5, 'rgb(254, 254, 189)'],
                            [0.6, 'rgb(217, 239, 139)'],
                            [0.7, 'rgb(164, 216, 105)'],
                            [0.8, 'rgb(102, 189, 99)'],
                            [0.9, 'rgb(63, 170, 89)'],              
                            [1.0, 'rgb(25, 151, 79)']]

        self.tree = Phylo.read(self.path, self.format)
        self.treeList = self.calculate(self.tree)
        #self.tree = Phylo.read(self.path, self.format)
        
        '''
        self.nleaves = self.tree.count_terminals()
        self.nnodes = len(self.tree.depths())
        self.nroots = self.nnodes - self.nleaves
        self.leaves = self.tree.get_terminals()
        self.roots = self.tree.get_nonterminals()
        '''

    def calculate(self, tree):
        traverse_order = 'preorder'
        all_clades = list(tree.find_clades(order=traverse_order))

        start_degree = self.degreerange[0] * np.pi/180
        end_degree = self.degreerange[1] * np.pi/180
        
        degreerange = 2*np.pi - (end_degree - start_degree)

        n = 0
        for k in range(len(all_clades)):
            all_clades[k].id = k
            
            if all_clades[k].is_terminal():
                n += 1
                all_clades[k].theta = start_degree + degreerange * n/self.tree.count_terminals()

        for k in reversed(range(len(all_clades))):
            if not all_clades[k].is_terminal():
                all_clades[k].theta = (all_clades[k][0].theta + all_clades[k][-1].theta)/2.0


        for k in range(len(all_clades)):

            # no need for .color, as we'll make it in the trace(marker=(colorbar=...))

            all_clades[k].radius = self.get_branch_length_dict()[all_clades[k]]

            all_clades[k].complex = np.zeros(1, dtype='complex')
            all_clades[k].complex.real = np.sin(all_clades[k].theta) * all_clades[k].radius
            all_clades[k].complex.imag = np.cos(all_clades[k].theta) * all_clades[k].radius

            if all_clades[k].confidence:
                all_clades[k].hovertext = 'id: {}<br>branch-length: {}<br>confidence: {}'.format(all_clades[k].id, all_clades[k].branch_length, all_clades[k].confidence.value)
                
            else:
                all_clades[k].hovertext = 'id: {}<br>name: {}<br>branch-length: {}'.format(all_clades[k].id, all_clades[k].name, all_clades[k].branch_length)
                
        return all_clades

    def complexList(self, radius=True):

        complexList = []
        if radius:
            for k in range(len(self.treeList)):
                complexList.append(np.zeros(1, dtype='complex'))
                complexList[k].real = np.sin(self.treeList[k].theta) * self.treeList[k].radius
                complexList[k].imag = np.cos(self.treeList[k].theta) * self.treeList[k].radius
        else:
            for k in range(len(self.treeList)):
                complexList.append(np.zeros(1, dtype='complex'))
                try:
                    complexList[k].real = np.sin(self.treeList[k].theta) * (self.treeList[k].radius - self.treeList[k].branch_length)
                    complexList[k].imag = np.cos(self.treeList[k].theta) * (self.treeList[k].radius - self.treeList[k].branch_length)
                except TypeError:
                    complexList[k].real = 0.0
                    complexList[k].imag = 0.0
        return complexList


    def color_size(self):
        color_list = []
        all_clades = self.calculate(self.tree)
        for k in range(len(all_clades)):
            if all_clades[k].confidence:
                color_list.append(all_clades[k].confidence.value)
            else:
                color_list.append(0.0)
        size_list = [9 if c!=-1 else 7 for c in color_list]
        
        return color_list, size_list

    def check_unit_branch_lengths(self):
        
        if self.cladogram is True:
            return True
        elif not np.count_nonzero(self.tree.depths().values()):
            return True
        else:
            return False

    def get_branch_length_dict(self):
        
        return self.tree.depths(unit_branch_lengths=self.check_unit_branch_lengths())


    def phylo_trace(self):
        
        x = [*map(lambda x:x.real.tolist(), self.complexList())]
        x = sum(x, [])
        y = [*map(lambda x:x.imag.tolist(), self.complexList())]
        y = sum(y, [])

        color, size = self.color_size()

        trace = [go.Scatter(x=x, 
                           y=y,
                           mode='markers',
                           hoverinfo='text',
                           text=[*map(lambda x:x.hovertext, self.treeList)],
                           opacity=1,
                           marker=dict(color=color,
                                       size=size, 
                                       colorscale=self.colorscale,
                                       colorbar=dict(thickness=20, dtick=10, ticklen=4, title='confidence')
                                       ))]
        
        return trace
    
    def get_radial_line(self):
        x = [*map(lambda x:x.real.tolist(), self.complexList())]
        x = sum(x, [])
        y = [*map(lambda x:x.imag.tolist(), self.complexList())]
        y = sum(y, [])

        x0 = [*map(lambda x:x.real.tolist(), self.complexList(radius=False))]
        x0 = sum(x0, [])
        y0 = [*map(lambda x:x.imag.tolist(), self.complexList(radius=False))]
        y0 = sum(y0, [])

        length = len(x)
        assert len(x) == len(y) == len(x0) == len(y0)
        path_list = np.column_stack((np.full(length,'M'), x, y, np.full(length, 'L'), x0, y0))
        path = " ".join(path_list.flatten())
        return path
    
    
    '''
    def get_arc(self):
         for k in reversed(range(len(all_clades))):
            if not all_clades[k].is_terminal():
                # all_clades[k].theta = (all_clades[k][0].theta + all_clades[k][-1].theta)/2.0
    '''
    
    
    def phylo_layout(self):
        layout = go.Layout(autosize=False,
                        showlegend=False,
                        xaxis=dict(visible=False),
                        yaxis=dict(visible=False),
                        hovermode='closest',
                        width=1500,
                        height=1500)
        layout['shapes'] = [dict(path=self.get_radial_line())]
        
        return layout

    def phylo_fig(self):
        return go.Figure(data=self.phylo_trace(), layout=self.phylo_layout())
        
        '''
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
            

            
            for i in range(len(node_complex)):
                trace.append(go.Scatter(x=node_complex[i].real,
                                        y=node_complex[i].imag,
                                        marker=dict(color='red', size=10),
                                        opacity=1))
            
            
            return trace

        def phylo_layout(self):
            layout = go.Layout(autosize=False,
                            showlegend=False,
                            xaxis=dict(visible=False),
                            yaxis=dict(visible=False),
                            hovermode='closest',
                            width=1500,
                            height=1500)
            #layout['shapes'] = [dict(path=self.get_stem_path()), dict(path=self.get_angular_path())]
            
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

    print(radius)
    print('\t')
    print([*map(lambda t:len(t), theta)])
    print('\t')
    print(theta)
    print('\t')
    print(complex_range)
    print('\t')


    print(path)
    print('\t')
    print(phylo.phylo_trace())



    print('total run time:')
    print (time()-t)




    fig = phylo.phylo_fig()
    plot(fig)
    '''




phylo = circularTree(path=sys.argv[1])
#
'''
print('phylo.phylo_trace() is:')
print('\t')
print(phylo.phylo_trace())
print('\t')
print('phylo.tree is:')
print('\t')
print(phylo.tree)
print('\t')
print('phylo.treeList is:')
print('\t')
print(phylo.treeList)
print('\t')
print('phylo.complexList is:')
print('\t')
print(phylo.complexList())
'''
print(phylo.complexList())
print('\t')


plot(phylo.phylo_fig())
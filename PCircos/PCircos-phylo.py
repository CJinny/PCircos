import sys
import plotly.graph_objs as go
from plotly.offline import download_plotlyjs, init_notebook_mode, plot, iplot
init_notebook_mode(connected=True)
import numpy as np
import pandas as pd
from Complex import Complex
from Bio import Phylo
from time import time
import json
from fig import Figure
from config import json_config
import colorlover as cl

# Dependency:
## Biopython

# important function in biopython Phylo:
# tree.depths()

def json_dict(json_path):
    with open(json_path, 'r') as f:
        json_read = f.read()
    f.close()
    json_dict = json.loads(json_read)

class circularTree(Figure):
    def __init__(self, *args, **kwargs):
        ## TEMP
        assert 'config' in kwargs
        self.config = json_config.json_dict(kwargs.pop('config'))
        self.treepath = self.config['Category']['tree']['file']['path']
        self.treeformat = self.config['Category']['tree']['file']['format']
        self.cladogram = self.config['Category']['tree']['cladogram']
        self.confidence = self.config['Category']['tree']['confidence']
        self.degreerange = self.config['Category']['tree']['degreerange']
        self.npoints = self.config['Category']['tree']['npoints']
        self.start_leaf = self.config['Category']['tree']['start_leaf']
        if 'meta' in self.config['Category']['tree'] and self.config['Category']['tree']['meta']['meta']:
            self.meta = pd.read_csv(self.config['Category']['tree']['meta']['path'],
                                   header = self.config['Category']['tree']['meta']['header'],
                                   sep = self.config['Category']['tree']['meta']['sep'])
        else:
            self.meta = None

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

        self.tree = Phylo.read(self.treepath, self.treeformat)
        self.treeList = self.fill_tree(self.tree)
         

    def fill_tree(self, tree):
        # in cases where only confidence value color needs to be drawn, self.confidence2color

        traverse_order = 'preorder'
        all_clades = list(tree.find_clades(order=traverse_order))

        start_degree = self.degreerange[0] * np.pi/180
        end_degree = self.degreerange[1] * np.pi/180
        
        degreerange = 2*np.pi - (end_degree - start_degree)

        n = 0
        for k in range(len(all_clades)):
            
            all_clades[k].id = k

   

            if all_clades[k].is_terminal():
                if self.meta is not None:
                    all_clades[k].color = self.meta2color()[all_clades[k].name]
                n += 1
                all_clades[k].theta = start_degree + degreerange * n/self.tree.count_terminals()
              

        for k in reversed(range(len(all_clades))):
            if not all_clades[k].is_terminal():
                all_clades[k].theta = (all_clades[k][0].theta + all_clades[k][-1].theta)/2.0

                if self.meta is not None:
                    
                    leaves = [*map(lambda t: self.meta2color()[t.name], all_clades[k].get_terminals())]
                    print('leaves are:')
                    print(leaves)

                    all_clades[k].color = tuple(np.mean(np.array(leaves)), axis=0)
                    #all_clades[k].color = 

        for k in range(len(all_clades)):

            all_clades[k].radius = self.get_branch_length_dict()[all_clades[k]]
            # the first element in all_clades will be the center of the circle, sometimes the branch_length could be Nontype
            if k == 0:
                all_clades[k].root_radius = 0.0
            else:
                all_clades[k].root_radius = round(all_clades[k].radius - all_clades[k].branch_length, 6)

            all_clades[k].complex = np.zeros(1, dtype='complex')
            all_clades[k].complex.real = np.sin(all_clades[k].theta) * all_clades[k].radius
            all_clades[k].complex.imag = np.cos(all_clades[k].theta) * all_clades[k].radius

            if all_clades[k].confidence:
                all_clades[k].hovertext = 'id: {}<br>branch-length: {}<br>confidence: {}'.format(all_clades[k].id, round(all_clades[k].branch_length, 4), all_clades[k].confidence.value)
                
            else:
                all_clades[k].hovertext = 'id: {}<br>name: {}<br>branch-length: {}'.format(all_clades[k].id, all_clades[k].name, all_clades[k].branch_length)

        return all_clades

    def complex_array(self, radius=True):

        complex_array = np.zeros(len(self.treeList), dtype='complex')

        if radius:
            for k in range(len(self.treeList)):
                
                complex_array.real[k] = np.sin(self.treeList[k].theta) * self.treeList[k].radius
                complex_array.imag[k] = np.cos(self.treeList[k].theta) * self.treeList[k].radius
        else:
            for k in range(len(self.treeList)):
                try:
                    complex_array.real[k] = np.sin(self.treeList[k].theta) * (self.treeList[k].radius - self.treeList[k].branch_length)
                    complex_array.imag[k] = np.cos(self.treeList[k].theta) * (self.treeList[k].radius - self.treeList[k].branch_length)
                except TypeError:
                    complex_array.real[k] = 0.0
                    complex_array.imag[k] = 0.0
        return complex_array


    def confidence2col(self):
        if not self.confidence:
            return None
        else:
            color_list = []
            all_clades = self.treeList
            for k in range(len(all_clades)):
                if all_clades[k].confidence:
                    color_list.append(all_clades[k].confidence.value)
                else:
                    color_list.append(0.0)
            #size_list = [9 if c!=-1 else 7 for c in color_list]
            
            return color_list

    def meta2color(self):
        # ONGOING
        if self.meta is None:
            return None
        else:
            color_keys = self.meta.iloc[:,self.config['Category']['tree']['meta']['idcolumn']]
            color_categories = self.meta.iloc[:,self.config['Category']['tree']['meta']['colorcolumn']]
            
            unique_categories = np.unique(color_categories)

            warm = cl.scales['5']['seq']['YlOrBr']
            cold = cl.scales['5']['seq']['YlGnBu']
            length = round(len(unique_categories)/2)
            
            t = 5 - length % 5 
            length += t
            d = length * 2 - len(unique_categories)


            unique_color_vals = list(reversed(cl.interp(cold,length)))
            unique_color_vals.extend(cl.interp(warm,length))
            unique_color_vals = cl.to_numeric(cl.to_rgb(unique_color_vals))

            while len(unique_color_vals) > len(unique_categories):
                unique_color_vals.pop(len(unique_color_vals)//2)
   
            unique_color_dict = dict(zip(unique_categories, unique_color_vals))

            color_vals = [*map(lambda t:unique_color_dict[t], color_categories)]
            
            '''
            print('\t')
            print('color_vals is:\t')
            print(color_vals)
            print('\t')
            '''
            return dict(zip(color_keys,color_vals))


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
        
       
        x, y = self.complex_array().real, self.complex_array().imag

        if self.confidence:
            colors = self.confidence2col()
        
        elif self.meta is not None:
            colors = [*map(lambda t:t.color, self.treeList)]

        else:
            colors = np.repeat('black',len(x))

        trace = [go.Scatter(x=x, 
                            y=y,
                            mode='markers',
                            hoverinfo='text',
                            text=[*map(lambda x:x.hovertext, self.treeList)],
                            opacity=1,
                            marker=dict(color=colors,
                                        size=8, 
                                        #colorscale=self.colorscale,
                                        #colorbar=dict(thickness=20, dtick=10, ticklen=4, title='confidence')
                                        ))]
        
        return trace
    
    def get_radial_line(self):


        x, y = self.complex_array().real, self.complex_array().imag
        x0, y0 = self.complex_array(radius=False).real, self.complex_array(radius=False).imag

        arcid = []
        def generate_arcid(clade):
            tmp = []
            for i in range(len(clade)):
                tmp.append(clade[i].id)
            arcid.append(tmp)
            if not clade.is_terminal():
                for n in range(len(clade)):
                    if not clade[n].is_terminal():
                        generate_arcid(clade[n])
        generate_arcid(self.tree.clade)
  

        def map_id(id, list=self.treeList, attribute='theta'):
            res = []
            for i in range(len(id)):
                res.append([])
                if attribute == 'theta':
                    res[i] = [*map(lambda t: list[t].theta, id[i])]
                elif attribute == 'root_radius':
                    res[i] = [*map(lambda t: list[t].root_radius, id[i])]
            return res

        theta_list = map_id(arcid)

        #theta_list = np.array([*map(lambda x: x.theta, self.treeList)])[arcid]

        theta_linspace = [*map(lambda x: np.linspace(x[0],x[-1],int(self.npoints * (np.max(x)-np.min(x))/(2*np.pi))+2), theta_list)]

        arc_radius = map_id(arcid, attribute='root_radius')

        #arc_radius = np.array([*map(lambda x: x.root_radius, self.treeList)])[arcid]
       
        arc_complex = [*map(lambda x, y: np.sin(x)*np.mean(y) + np.array([1j])*np.cos(x)*np.mean(y),  theta_linspace, arc_radius)]

        
        # the circle center, we can pop that
        arc_complex.pop(0)
        arc_path_list = [*map(lambda x: " ".join(np.column_stack((np.concatenate((np.full(1,'M'), np.full(len(x)-1,'L'))), x.real, x.imag)).ravel()), arc_complex)]


        length = len(x)
        assert len(x) == len(y) == len(x0) == len(y0)
        radial_path_list = np.column_stack((np.full(length,'M'), x, y, np.full(length, 'L'), x0, y0)).ravel().tolist()
        
        path_list = arc_path_list + radial_path_list
        
        path = " ".join(path_list)
        return path
   
    
    
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
        
 

t=time()

phylo = circularTree(config=sys.argv[1])
'''
if len(sys.argv) -1 == 3:
    phylo = circularTree(path=sys.argv[1], format=sys.argv[2])

elif len(sys.argv) -1 == 2:
    phylo = circularTree(path=sys.argv[1], format=sys.argv[2])
else:
    phylo = circularTree(path=sys.argv[1])
'''


plot(phylo.phylo_fig())
print (time()-t)
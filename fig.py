
import sys
import plotly.graph_objs as go
import numpy as np
import pandas as pd
import maths
import colors
import copy
from Complex import Complex
from config import json_config, coord_config


# input_json_path = sys.argv[1]

## TO DO:
    # ABORT merge_dict()
    # ABORT customcolor in annotation, this is so weird!



def merge_dict(basedict, *extradict):
    '''
    this function updates basedict with extradict values, will apply to any nested values
    e.g:
    extradict = {'font': {'color': 'yellow'}}
    basedict = {'font': {'color': 'black', size: 16}}
    merge_dict(extradict, basedict):
        {'font': {'color': 'yellow', size: 16}}
    '''
    for i in range(len(extradict)):
        for key in extradict[i].keys():
            if key not in basedict.keys():
                basedict[key] = extradict[i][key]
            else:
                if not isinstance(basedict[key], dict):
                    basedict[key] = extradict[i][key]

                else:
                    merge_dict(basedict[key], extradict[i][key])
  
    return basedict

class Figure(Complex):
    def __init__(self, *args, **kwargs):
        
        assert 'input_json_path' in kwargs
        if 'dash_dict' in kwargs:
            self.dash_dict = kwargs.pop('dash_dict')
        else:
            self.dash_dict = None


        self.config_dict = json_config.json2dict(kwargs.pop('input_json_path'), dash_dict=self.dash_dict)




        #super(Figure, self).__init__(*args, **kwargs)

        self.ideogram_coord_config = self.get_ideogram_coord_config()

        self.layout_general = self.config_dict['General']

        self.categories = self.config_dict['Category']

        self.ideogram = self.categories['ideogram']

        self.ideogram_ideogram = self.ideogram['ideogram']


        self.ideogram_majortick = self.ideogram['majortick']
        self.ideogram_minortick = self.ideogram['minortick']
        self.ideogram_ticklabel = self.ideogram['ticklabel']


        self.ideogram_radius_dict = self.ideogram_ideogram['radius']
        self.show_chr_annotation = self.ideogram_ideogram['chrannotation']['show']
        #self.show_major_tick = self.ideogram_majortick['show']
        #self.show_minor_tick = self.ideogram_minortick['show']

        
        self.major_tick_radius_dict = self.ideogram_majortick['radius']
        self.minor_tick_radius_dict = self.ideogram_minortick['radius']
        self.tick_label_radius_dict = self.ideogram_ticklabel['radius']
        self.SUM = self.get_ideogram_coord_config()['SUM']
        self.degreerange = self.ideogram_ideogram['degreerange']
        self.chr_color_dict = self.get_chr_info()['chr_color_dict']   # a dict: chr_name:chr_color
        self.chr_label_dict = self.get_chr_info()['chr_label_dict']


    def np_list_concat(self, x):
        if isinstance(x, list):
            try:
                return np.concatenate(x)
            except ValueError:
                return np.array(x)

        elif isinstance(x, np.ndarray):
            return x
        else:
            raise ValueError('input must be an ndarray or a list')

    def get_read_data(self, key):
        
        items = self.categories[key]
        sortindices = self.get_data_array_sortindex(key)


        def get_single_data(item, key, sortindex=None):
            
            assert isinstance(item, dict)



            if item['file']['header'] in ['None', None, 'none']:
                unsorted_data = coord_config.read_data(item['file']['path'],
                                                       key,
                                                       self.get_chr_info(),
                                                       sep=item['file']['sep'],
                                                       header=None
                                                       )
            else:
                unsorted_data =  coord_config.read_data(item['file']['path'],
                                                        key,
                                                        self.get_chr_info(),
                                                        sep=item['file']['sep']
                                                        )
            if item['sortbycolor']:
                ## ONGOING
                assert sortindex is not None
                return unsorted_data[sortindex]
            else:
                return unsorted_data
            
        if isinstance(items, dict):
            return get_single_data(items, key, sortindex=sortindices)
        elif isinstance(items, list):
            return [*map(lambda x, y: get_single_data(x, key, sortindex=y), items, sortindices)]   


    def get_chr_info(self):
        # temp
        #chr_info_file = self.ideogram_ideogram['file']
        #custom_options = self.ideogram_ideogram['customoptions']

        chr_info_file = self.config_dict['Category']['ideogram']['ideogram']['file']
        custom_options = self.config_dict['Category']['ideogram']['ideogram']['customoptions']

        if chr_info_file['header'] in ['None', None, 'none']:
            chr_info_dict = coord_config.chr_info(chr_info_file['path'], 
                                                  sep=chr_info_file['sep'], 
                                                  header=None,
                                                  custom_label=custom_options['customlabel'], 
                                                  custom_spacing=custom_options['customspacing'], 
                                                  custom_color=custom_options['customcolor'],
                                                  dash_dict=self.dash_dict
                                                  )
        else:
            chr_info_dict = coord_config.chr_info(chr_info_file['path'], 
                                                  sep=chr_info_file['sep'], 
                                                  custom_label=custom_options['customlabel'], 
                                                  custom_spacing=custom_options['customspacing'], 
                                                  custom_color=custom_options['customcolor'],
                                                  dash_dict=self.dash_dict
                                                  )
        return chr_info_dict                                          

        
    def get_ideogram_coord_config(self):
        # temp comment out

        ideogram_coord_config = coord_config.ideogram_coord_config(self.get_chr_info(), 
                                                                   #npoints=self.ideogram_ideogram['npoints'],
                                                                   npoints=self.config_dict['Category']['ideogram']['ideogram']['npoints'],
                                                                   show_major_tick=self.config_dict['Category']['ideogram']['majortick']['show'],
                                                                   major_tick_spacing=self.config_dict['Category']['ideogram']['majortick']['spacing'],
                                                                   show_minor_tick=self.config_dict['Category']['ideogram']['minortick']['show'],
                                                                   minor_tick_spacing=self.config_dict['Category']['ideogram']['minortick']['spacing'],
                                                                   show_tick_label=self.config_dict['Category']['ideogram']['ticklabel']['show'],
                                                                   tick_label_spacing=self.config_dict['Category']['ideogram']['ticklabel']['spacing']
                                                                  )
        return ideogram_coord_config

    ## self.ideogram_complex() and self.ring_complex() is inherited
   
    def get_ideogram_theta(self):
        return self.ideogram_theta_list(self.ideogram_coord_config, self.SUM, degreerange=self.degreerange)

    def get_ideogram_complex(self):
        return self.ideogram_complex(self.ideogram_coord_config, self.SUM, self.degreerange, ideogram_radius_dict=self.ideogram_radius_dict)

    #def get_ideogram_shapes(self):
        #return self.ideogram_path(self.get_ideogram_complex())
       

    def get_ring_complex(self):
        if isinstance(self.categories['ring'], dict):
            return self.ideogram_complex(self.ideogram_coord_config, self.SUM, degreerange=self.degreerange, ideogram_radius_dict=self.categories['ring']['radius'])
        elif isinstance(self.categories['ring'], list):
             return [*map(lambda x: self.ideogram_complex(self.ideogram_coord_config, self.SUM, degreerange=self.degreerange, ideogram_radius_dict=x['radius']), self.categories['ring'])]

    def get_ring_paths_dict(self):
        
        def single_ring_dict(path, ring_dict):
            # ring data can only have one background color! therefore I join them
            # item = self.categories['ring']
            if isinstance(path, list):
                path = self.pathjoin(path)
            path_dict=dict(path=path)
            path_dict.update(ring_dict['layout'])
            return path_dict

        if isinstance(self.categories['ring'], dict):
            path = self.ideogram_path(self.get_ring_complex())
            path_dict =[single_ring_dict(path, self.categories['ring'])]

        elif isinstance(self.categories['ring'], list):
            path = [*map(lambda x: self.ideogram_path(x), self.get_ring_complex())]
            path_dict = [*map(lambda x, y: single_ring_dict(x, y), path, self.categories['ring'])]
        
        return path_dict


    def get_major_tick_path(self):
        major_tick_accum_coord_list = self.ideogram_coord_config['major_tick_accum_coord_list']
        major_tick_theta = self.ideogram_tick_theta_list(self.ideogram_coord_config, major_tick_accum_coord_list, SUM=self.SUM, degreerange=self.degreerange)
        major_tick_complex = self.tick_complex(major_tick_theta, tick_radius_dict=self.major_tick_radius_dict)

        return self.tick_path(major_tick_complex)
       
    def get_minor_tick_path(self):
        minor_tick_accum_coord_list = self.ideogram_coord_config['minor_tick_accum_coord_list']
        minor_tick_theta = self.ideogram_tick_theta_list(self.ideogram_coord_config, minor_tick_accum_coord_list, SUM=self.SUM, degreerange=self.degreerange)
        minor_tick_complex = self.tick_complex(minor_tick_theta, tick_radius_dict=self.minor_tick_radius_dict)

        return self.tick_path(minor_tick_complex)   

    def get_ideogram_chrannot_theta(self):
        return self.ideogram_chrannot_theta(self.ideogram_coord_config, self.SUM, degreerange=self.degreerange)

    def get_ideogram_chrannot_complex(self):
        return self.ideogram_chrannot_complex(self.SUM,
                                              degreerange=self.degreerange,
                                              chr_annotation_radius_dict=self.ideogram_ideogram['chrannotation']['radius'])
    
    def get_tick_label_complex(self):
        tick_label_theta = self.ideogram_tick_label_theta_list(self.ideogram_coord_config, self.SUM, degreerange=self.degreerange)
        return self.tick_label_complex(tick_label_theta, tick_label_radius_dict=self.tick_label_radius_dict)
        

    def pathjoin(self, path_list):
        '''this function will join any path_string_list into a single path_string, useful when fillcolor is the same'''
        return " ".join(path_list)

    def get_data_array_dict(self, key):
        # single instance key
        # without ideogram
        # for each category type
        ## either returns a np.ndarray or a list of np.ndarray!
        assert key in self.categories
        items = self.categories[key]
        if key == 'ideogram':
            raise ValueError('ideogram information should not be parsed in get_data_array()')
        else:
 
            def single_data_array(key, item):

                if 'colorcolumn' not in item.keys():
                    item['colorcolumn'] = None  
                if 'sortbycolor' not in item.keys():
                    item['sortbycolor'] = False

                if item['file']['header'] in ['None', None, 'none']:
                    return coord_config.data_array(item['file']['path'], key,
                                                self.get_chr_info(), 
                                                sep=item['file']['sep'],
                                                header=None,
                                                colorcolumn=item['colorcolumn'],
                                                sortbycolor=item['sortbycolor']
                                                )
                else:
                    return coord_config.data_array(item['file']['path'], key,
                                                self.get_chr_info(), 
                                                sep=item['file']['sep'],
                                                colorcolumn=item['colorcolumn'],
                                                sortbycolor=item['sortbycolor']
                                                )

            if isinstance(items, dict):
                return single_data_array(key, items)
            elif isinstance(items, list):
                return [*map(lambda x: single_data_array(key, x), items)]

    def get_data_array(self, key):
        if not isinstance(self.get_data_array_dict(key), list):
            return self.get_data_array_dict(key)['data_array']
        else:
            return [*map(lambda x: x['data_array'], self.get_data_array_dict(key))]

    def get_data_array_sortindex(self, key):
        # only used when sortbycolor is True! in other cases its just range(len(data_array))

        if not isinstance(self.get_data_array_dict(key), list):
            return self.get_data_array_dict(key)['sortindex']
        else:
            return [*map(lambda x: x['sortindex'], self.get_data_array_dict(key))]


    def get_data_complexes(self, key, return_path=True):
        assert key in self.categories
        data_array = self.get_data_array(key)
        items = self.categories[key]

        def single_data_complex(data_array, key, item):
            

            if key != 'highlight':
                
                if key == 'cytoband':
                    item['radius'] = self.ideogram_radius_dict
                elif key == 'annotation' and item['customradius']:
                    try:
                        assert isinstance(item['radiuscolumn'], int)
                    except AssertionError:
                        print ('Please enter a valid radiuscolumn under annotation')
                    try:
                        radiuscolumn = item['radiuscolumn']
                    except IndexError:
                        print ('the column you entered is out of bound, notice that column starts with 0, not 1 ')
                    try:
                        data_array[:,radiuscolumn].astype('float')
                    
                    ## DEBUG, the below should be ValueError, changing to Exception temporarily
                    except Exception:
                        print (radiuscolumn)
                        print (data_array[:,radiuscolumn])
                        print ('Please make sure to enter numeric value for radius column')
                
                    item['radius'] = {"R": data_array[:,radiuscolumn]}
                
                radius_dict = item['radius']
                
            else:

                radius_dict = {"R0": data_array[:,item['R0column']], "R1": data_array[:,item['R1column']]}
                

            if key == 'annotation':
                data_complex = self.data_complex(self.ideogram_coord_config, 
                                                 data_array, key, radius_dict, self.SUM, 
                                                 degreerange=self.degreerange,
                                                 custom_offset_degree=item['customoffsetdegree']
                                                 )
            else:
                
                data_complex = self.data_complex(self.ideogram_coord_config, 
                                                 data_array, key, radius_dict, self.SUM, 
                                                 degreerange=self.degreerange,
                                                 return_path=return_path
                                                 )
            return data_complex
        if isinstance(items, dict):
            return single_data_complex(data_array, key, items)
        elif isinstance(items, list):
            return [*map(lambda x, y: single_data_complex(x, key, y), data_array, items)]




    def get_hovertext(self, key):

        assert key in ['histogram', 'line', 'area','scatter', 'tile', 'heatmap', 'link', 'ribbon', 'twistedribbon']

        if not isinstance(self.categories[key], list):
            hovertextformat = self.categories[key]['hovertextformat']

            a = self.get_read_data(key)
            hvtext = []
            if key not in ['link', 'ribbon', 'twistedribbon']:
                assert a.shape[1] >= 3

                if key in ['histogram', 'tile', 'heatmap']:
                    
                    list_count = [*map(lambda x: len(x), self.get_data_complexes(key))]

                    for i in range(len(a)):
                        k = 0
                        while k < list_count[i]:
                            k += 1
                            hvtext.append(eval(hovertextformat))                          
                else:
                    for i in range(len(a)):
                        hvtext.append(eval(hovertextformat))
                        
            else:
                # for link, ribbon and twistedribbon hovertext is a list of two elememnt
                assert a.shape[1] >= 6
                for i in range(len(a)):
                    hvtext.append(eval(hovertextformat[0]))
                    hvtext.append(eval(hovertextformat[0]))
                    hvtext.append(eval(hovertextformat[1]))
                    hvtext.append(eval(hovertextformat[1]))
        else:
            hvtext = []
            data_array_list = self.get_read_data(key)
            for j in range(len(data_array_list)):
                hovertextformat = self.categories[key][j]['hovertextformat']
                
                hvtext.append([])
                a = data_array_list[j]

                if key not in ['link', 'ribbon', 'twistedribbon']:
                    assert a.shape[1] >= 3
                    
                    if key in ['histogram', 'tile', 'heatmap']:
                        
                        list_count = [*map(lambda x: len(x), self.get_data_complexes(key)[j])]
                        
                        for i in range(len(a)):
                            k = 0
                            while k < list_count[i]:
                                k += 1
                                hvtext[j].append(eval(hovertextformat))

                    else: 
                        for i in range(len(a)):
                            hvtext[j].append(eval(hovertextformat))

                else:
                    assert a.shape[1] >= 6
                    for i in range(len(a)):
                        hvtext[j].append(eval(hovertextformat[0]))
                        hvtext[j].append(eval(hovertextformat[0]))
                        hvtext[j].append(eval(hovertextformat[1]))
                        hvtext[j].append(eval(hovertextformat[1]))
        return hvtext

    def get_traces(self, key):
        ## always return a list, this makes concatenation easier
        # for line plot, there shouldn't be sortbycolor

        items = self.categories[key]
        complexes = self.get_data_complexes(key, return_path=False)
        data_arrays = self.get_data_array(key)
        hovertexts = self.get_hovertext(key)


        def single_trace(key, Complex, item, data_array, hovertext):
            assert key not in ['cytoband', 'ideogram', 'ring', 'annotation', 'highlight', 'connector']
            # For scatter plot the complex is an ndarray, for lines it would be a list of ndarray separated by chromosomes
            # for other nonvisible plots, they can be concatenated into one ndarray
            assert isinstance(item, dict)

            if key != 'line':
                if isinstance(Complex, list):
                    Complex = np.concatenate(Complex)


            if key == 'line':
                # the only time when Complex is a list of ndarray
                trace = []
                index = np.cumsum([0]+[*map(lambda x: len(x), Complex)])

                def divide(l, index):
                    '''this function creates a generator object for hovertext so I know how many hovertext element to take for each chromosome'''
                    for n in range(len(index)-1):
                        yield l[index[n]:index[n+1]]

                hovertext_generator = divide(hovertext, index)

                for i in range(len(Complex)):
                    trace.append(go.Scatter(x=Complex[i].real,
                                            y=Complex[i].imag,
                                            text=next(hovertext_generator)
                                            )
                                )
                    trace[i].update(item['trace'])


            else:
                if Complex.ndim != 1:
                    Complex = Complex.ravel()
                    
                if key in ['link', 'ribbon', 'twistedribbon']:
                    index = []
                    for i in range(len(Complex)):
                        if i%6 in [0,1,4,5]:
                            index.append(i)
                    Complex = Complex[index]

                # all cases except when key == 'line'
                trace = go.Scatter(x=Complex.real,
                                   y=Complex.imag,
                                   text=hovertext
                                   )
                trace.update(item['trace'])

                if key == 'scatter':
                    if 'color' not in item['trace']['marker']:

                        chr_label = data_array[:,0]

                        if item['colorcolumn'] is None:
                            

                            color = [*map(lambda x: self.chr_color_dict[x], chr_label)]
                        else:
                            n = item['colorcolumn']
                            assert isinstance(n, int)
                            color = colors.to_rgb(data_array[:,n])   
                    else:
                        color = item['trace']['marker']['color']

                    trace['marker'].update(color=color)

            return trace

        if isinstance(items, dict):
            trace = single_trace(key, complexes, items, data_arrays, hovertexts)
            if not isinstance(trace, list):
                trace = [trace]
            return trace

        elif isinstance(items, list):
            return [*map(lambda w, x, y, z: single_trace(key, w, x, y, z), complexes, items, data_arrays, hovertexts)]



    def trace(self):
        # aggregate all get_traces element into one trace, list variable
        trace = []
        for key in self.categories.keys():
            if key not in ['cytoband', 'ideogram', 'ring', 'annotation', 'highlight', 'connector']:

                trace += self.get_traces(key)

        return trace

    def get_paths_dict(self, key):
        # will join path_list into a path string if sortbycolor
        ### deal with: 
            # histogram, ribbon, twistedribbon, if fillcolor is true, then linecolor==fillcolor!
            # cytoband, heatmap (sortbycolor=True)
            # tile, link (no fillcolor), color indicates to linecolor
            # area must have the same background area color! sortbycolor is disabled

        ## 
        assert key not in ['scatter', 'annotation', 'line']
        items = self.categories[key]
        data_arrays = self.get_data_array(key)
        data_complexes = self.get_data_complexes(key)

        def single_path(key, data_array, data_complex, item):
           
            if key in ['ribbon', 'twistedribbon']:
                interval_theta_array_0 = maths.to_theta(data_array[:,1:3], self.SUM, degreerange=self.degreerange)
                interval_theta_array_1 = maths.to_theta(data_array[:,4:6], self.SUM, degreerange=self.degreerange)
            else:
                interval_theta_array_0, interval_theta_array_1 = None, None

            if key == 'highlight':
                path_list = self.data_path(self.ideogram_coord_config, 
                                           key, data_complex, self.SUM, degreerange=self.degreerange,
                                           radius_dict={"R0": item['R0column'], "R1": item['R1column']},
                                           interval_theta_array_0=interval_theta_array_0,
                                           interval_theta_array_1=interval_theta_array_1)
            else:
                if key == 'cytoband':
                    item['radius'] = self.ideogram_radius_dict
                path_list = self.data_path(self.ideogram_coord_config, 
                                           key, data_complex, self.SUM, degreerange=self.degreerange,
                                           radius_dict=item['radius'],
                                           interval_theta_array_0=interval_theta_array_0,
                                           interval_theta_array_1=interval_theta_array_1)


            if not (item['sortbycolor'] and item['colorcolumn']):
                path = " ".join(path_list)
                paths_dict = dict(path=path)
                paths_dict.update(item['layout'])
            else:
                # custom line color only, and no fill color: link, tile
                # custom fill color only: all others

                paths_dict = []
                n = item['colorcolumn']
                color = data_array[:,n]

                if key == 'heatmap':

                    color = maths.val2heatmap(color, palatte_dict=item['palatte'])


                if key == 'highlight':
                    o = item['opacitycolumn']
                    opacity = data_array[:,o]


                    for i in range(len(color)):
                        paths_dict.append(dict(path=path_list[i],
                                               fillcolor=color[i],
                                               opacity=opacity[i]
                                               )
                                        )
                        
                        paths_dict[i].update(item['layout'])
                    

                elif key in ['heatmap', 'cytoband', 'ribbon', 'twistedribbon', 'tile', 'link', 'area']:

                    for i in range(len(path_list)):
                        paths_dict.append(dict(path=path_list[i]))
                        ## ONGOING, copy.deepcopy(item['layout'])
                        paths_dict[i].update(copy.deepcopy(item['layout']))
                        ##########################
                        #  ONGOING, testing heatmap
                        ## problem with importing from json
                        if key in ['tile', 'link']:
                            paths_dict[i]['line']['color'] = color[i]
                        else:
                            paths_dict[i]['line']['color'] = color[i]
                            paths_dict[i]['fillcolor'] = color[i]

            return paths_dict

        if isinstance(items, dict):
            if isinstance(single_path(key, data_arrays, data_complexes, items), dict):
                return [single_path(key, data_arrays, data_complexes, items)]
            else:
                return single_path(key, data_arrays, data_complexes, items)
            
        elif isinstance(items, list):
            ## ONGOING
            return [*map(lambda x, y, z: single_path(key, x, y, z), data_arrays, data_complexes, items)]
   
        else:
            raise KeyError('{} file does not exist!'.format(key))
    
   

    def get_annotations_dict(self):
        '''layout['annotations'] can only append one text annotation at a time!! be careful'''

        # please make sure you only have one annotation file!

        assert 'annotation' in self.categories.keys()
        try:
            assert isinstance(self.categories['annotation'], dict)
        except AssertionError:
            print ('Please use one and only one annotation file')

        text_complex = self.get_data_complexes('annotation')
        text_array = self.get_data_array('annotation')[:,2]
        text_theta = maths.to_theta(self.get_data_array('annotation')[:,1], self.SUM, degreerange=self.degreerange)

        
        textangle=self.angleconvert(text_theta, angleoffset=self.categories['annotation']['textangle']['angleoffset'],
                                    anglelimit=self.categories['annotation']['textangle']['anglelimit'])

        if self.categories['annotation']['customoffsetdegree']:
            textangle += self.get_data_array('annotation')[:,3]


        if self.categories['annotation']['fonttype'] == 'bold':
            if len(text_array) == 1:
                text = '<b>{}</b>'.format(text_array)
            else:
                text = [*map(lambda x: '<b>{}</b>'.format(x), text_array)]
        elif self.categories['annotation']['fonttype'] == 'italic':
            if len(text_array) == 1:
                text = '<i>{}</i>'.format(text_array)
            else:
                text = [*map(lambda x: '<i>{}</i>'.format(x), text_array)]
        elif self.categories['annotation']['fonttype'] in ['bold+italic', 'italic+bold']:
            if len(text_array) == 1:
                text = '<b><i>{}</i></b>'.format(text_array)
            else:
                text = [*map(lambda x: '<b><i>{}</i></b>'.format(x), text_array)]
        else:
            if len(text_array) == 1:
                text = '{}'.format(text_array)
            else:
                text = [*map(lambda x: '{}'.format(x), text_array)] 



        if len(text_array) == 1:
            
            annotation_dict = dict(x=text_complex.real, y=text_complex.imag,
                                   text=text, textangle=textangle)
            annotation_dict.update(copy.deepcopy(self.categories['annotation']['layout']))
            return [annotation_dict]
        else:
            if not self.categories['annotation']['customcolor']:
                annotation_dict_list =  [*map(lambda a,b,c: merge_dict(dict(x=a.real, y=a.imag, text=b, textangle=c), 
                                                                       copy.deepcopy(self.categories['annotation']['layout'])),
                                                            text_complex, text, textangle)]
            else:
                assert isinstance(self.categories['annotation']['colorcolumn'], int)
                n = self.categories['annotation']['colorcolumn']
                assert n <= self.get_data_array('annotation').shape[1]-1
                colors = self.get_data_array('annotation')[:,n].tolist()
                annotation_dict_list = [*map(lambda a,b,c,d: merge_dict(dict(x=a.real, y=a.imag, text=b, textangle=c), 
                                                                        copy.deepcopy(self.categories['annotation']['layout']),
                                                                        dict(font=dict(color=d))
                                                                        ),
                                                           text_complex, text, textangle, colors)]
            return annotation_dict_list
 


    def angleconvert(self, theta, angleoffset=-90, anglelimit=360, custom_offset_degree=0):
        # theta could be an ndarray or a list of ndarray
        # returns a processed degree data based on anglelimit
        # anglelimit is the degree limit before applying angleoffset and custom_offset!
        # custom_offset_degree is a 1D ndarray of the same shape as theta


        constant = 180/np.pi
        assert isinstance(angleoffset, (int, float))
        assert isinstance(anglelimit, (int, float))

        if not isinstance(theta, list):
            
            degree = theta*constant
          
            for i in range(len(degree)):    
                while degree[i] >= anglelimit:
                    degree[i] -= anglelimit

            degree += angleoffset + custom_offset_degree
            return degree
        else:
            assert isinstance(custom_offset_degree, list)
            assert len(theta) == len(custom_offset_degree)

            degree = [*map(lambda x: x*constant), theta]
            for i in range(len(degree)):
                for j in range(len(degree[i])):
                    while degree[i][j] >= anglelimit:
                        degree[i][j] -= anglelimit
                degree[i] += angleoffset + custom_offset_degree[i]
            return degree


    def layout(self):
        ## new version of plotly says layout['shapes'], ['annotations'] is a tuple, can't be used with list append
        ## I'm trying to add a temporary layout_shapes list variable collecting everything before assigning to layout['shapes']

        layout = go.Layout(self.layout_general)

        ## ONGOING
        layout_shapes = []
        layout_annotations = []

        if self.ideogram_ideogram['show']:
            
            if not self.ideogram_ideogram['showfillcolor']:
                # if no show fillcolor, there is no need to separate pathstring into list, hense the join
                # seems like for dictionary update, the dict variable needs to be defined first

                ideogram_pathstring = self.pathjoin(self.ideogram_path(self.get_ideogram_complex()))


                layout_dict = self.ideogram_ideogram['layout']
                layout_dict.update(dict(path=ideogram_pathstring))
                layout_shapes.append(layout_dict)
            else:
                ideogram_path_list = self.ideogram_path(self.get_ideogram_complex())

                for i in range(len(ideogram_path_list)):
                    layout_shapes.append(dict(path=ideogram_path_list[i], fillcolor=self.get_chr_info()['chr_fillcolor'][i]))
                    layout_shapes[i].update(self.ideogram_ideogram['layout'])
                    

        
        if self.ideogram_ideogram['chrannotation']['show']:
            
            chrannot_theta = self.get_ideogram_chrannot_theta()

            chrannot_theta = self.np_list_concat(chrannot_theta)

                
            chrannot_complex = maths.to_complex(chrannot_theta, self.ideogram_ideogram['chrannotation']['radius']['R'])

            chrannot_complex = self.np_list_concat(chrannot_complex)

            chrannot_angleoffset = self.ideogram_ideogram['chrannotation']['textangle']['angleoffset']
            chrannot_anglelimit = self.ideogram_ideogram['chrannotation']['textangle']['anglelimit']
                
            if self.ideogram_ideogram['chrannotation']['fonttype'] == 'bold':
                chrannot_text = [*map(lambda x: '<b>{}</b>'.format(x), self.get_chr_info()['chr_label'])]
            elif self.ideogram_ideogram['chrannotation']['fonttype'] == 'italic':
                chrannot_text = [*map(lambda x: '<i>{}</i>'.format(x), self.get_chr_info()['chr_label'])]
            elif self.ideogram_ideogram['chrannotation']['fonttype'] in ['bold+italic', 'italic+bold']:
                chrannot_text = [*map(lambda x: '<b><i>{}</i></b>'.format(x), self.get_chr_info()['chr_label'])]
            else:
                chrannot_text = '{}'.format(self.get_chr_info()['chr_label'])

            textangle = self.angleconvert(chrannot_theta, angleoffset=chrannot_angleoffset, anglelimit=chrannot_anglelimit)



            for i in range(len(chrannot_complex)):
                
                layout_annotations.append(dict(x=chrannot_complex[i].real,
                                                  y=chrannot_complex[i].imag,
                                                  text=chrannot_text[i],
                                                  textangle=textangle[i]
                                                  )
                                            )
                layout_annotations[i].update(self.ideogram_ideogram['chrannotation']['layout'])


        if self.ideogram_majortick['show']:
            
            layout_shapes.append(dict(path=self.get_major_tick_path()))
            layout_shapes[-1].update(self.ideogram_majortick['layout'])

        if self.ideogram_minortick['show']:
            
            layout_shapes.append(dict(path=self.get_minor_tick_path()))
            layout_shapes[-1].update(self.ideogram_minortick['layout'])

        if self.ideogram_ticklabel['show']:

            ticklabel_text = self.get_ideogram_coord_config()['tick_label_non_accum_list']
            if isinstance(ticklabel_text, list):
                ticklabel_text = np.concatenate(ticklabel_text)

            ticklabel_coord = self.np_list_concat(self.get_ideogram_coord_config()['tick_label_accum_coord_list'])
            ticklabel_theta = maths.to_theta(ticklabel_coord, self.SUM, degreerange=self.degreerange)
            ticklabel_angle = self.angleconvert(ticklabel_theta, angleoffset=self.ideogram_ticklabel['textangle']['angleoffset'], anglelimit=self.ideogram_ticklabel['textangle']['anglelimit'])
            
            ticklabel_complex = np.concatenate(self.get_tick_label_complex())

            if self.ideogram_ticklabel['textformat'] == 'Kb':
                ticklabel_text = [*map(lambda x: '{} Kb'.format(x//1000), ticklabel_text)]

            elif self.ideogram_ticklabel['textformat'] == 'Mb':
                ticklabel_text = [*map(lambda x: '{} Mb'.format(x//1000000), ticklabel_text)]

            elif self.ideogram_ticklabel['textformat'] == 'Gb':
                ticklabel_text = [*map(lambda x: '{} Gb'.format(x//1000000000), ticklabel_text)]

            else:
                raise ValueError('acceptable ideogram ticklabel textformats are: Kb & Mb & Gb')



            for i in range(len(ticklabel_complex)):
                layout_annotations.append(dict(x=ticklabel_complex[i].real,
                                                  y=ticklabel_complex[i].imag,
                                                  text=ticklabel_text[i],
                                                  textangle=ticklabel_angle[i]
                                                  )
                                            )
                layout_annotations[-1].update(self.ideogram_ticklabel['layout'])
            

   
            # due to the uncertain number of plot each type, we'll use extend [] instead of append!
        if 'ring' in self.categories.keys():
            # the idea is to always to draw the ring background first, whereas highlight and custom annotation will be drawn last!
               # ONGOING, maybe to deprecate show=True for ring?

            layout_shapes.extend(self.get_ring_paths_dict())



        for key in self.categories.keys():
            if key in ['histogram', 'cytoband', 'area', 'tile', 'heatmap', 'link','ribbon', 'twistedribbon', 'connector']:
                # ONGOING, maybe to deprecate self.categories[key]['show']
                #if self.categories[key]['show']:
                
                if isinstance(self.get_paths_dict(key)[0], dict):
                    try:
                        layout_shapes.extend(self.get_paths_dict(key))
                    except Exception:
                        print ('Error trying to plot {}'.format(key))
                        break
                elif isinstance(self.get_paths_dict(key)[0], list):
                    try:
                        layout_shapes.extend(sum(self.get_paths_dict(key),[]))
                    except Exception:
                        print ('Error trying to plot {}'.format(key))
                        break


        if 'highlight' in self.categories.keys():
            if self.categories['highlight']['show']:
                layout_shapes.extend(self.get_paths_dict('highlight'))

        if 'annotation' in self.categories.keys():
            if self.categories['annotation']['show']:
                layout_annotations.extend(self.get_annotations_dict())
        

        layout['shapes'] = layout_shapes
        layout['annotations'] = layout_annotations
        
        return layout



    def fig(self):
       
        return go.Figure(data=self.trace(), layout=self.layout())

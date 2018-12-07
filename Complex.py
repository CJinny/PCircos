import numpy as np
import maths, colors
import copy


# This module contains all information about Complex and path

class Complex(object):
    

        
    def __init__(self):
        pass

    def ideogram_tick_label_accum_coord_list(self, ideogram_coord_config):
        return ideogram_coord_config['tick_label_accum_coord_list']

    def ideogram_theta_list(self, ideogram_coord_config, SUM, degreerange=[0,360]):
        # 2X 
        return maths.to_theta(ideogram_coord_config['ideogram_accum_coord_list'], SUM, degreerange=degreerange)


    def ideogram_tick_theta_list(self, ideogram_coord_config, tick_accum_coord_list, SUM, degreerange=[0,360]):
        return maths.to_theta(tick_accum_coord_list, SUM, degreerange=degreerange)

    def ideogram_tick_label_theta_list(self, ideogram_coord_config, SUM, degreerange=[0,360]):
        return maths.to_theta(ideogram_coord_config['tick_label_accum_coord_list'], SUM, degreerange=degreerange)

    def ideogram_chrannot_theta(self, ideogram_coord_config, SUM, degreerange=[0,360]):
        return [*map(lambda x: np.mean(x), self.ideogram_theta_list(ideogram_coord_config, SUM, degreerange=degreerange))]
    
    def ideogram_chrannot_complex(self, SUM, degreerange=[0,360], chr_annotation_radius_dict={'R': 1.2}):
        return maths.to_complex(self.ideogram_chrannot_theta(SUM, degreerange), chr_annotation_radius_dict['R'])


    def tick_complex(self, tick_theta, tick_radius_dict={'R0': 1.1, 'R1': 1.12}):
        # used by both the major and minor tick
        return maths.to_complex(tick_theta, np.array([[tick_radius_dict['R0'], tick_radius_dict['R1']]]))

    def tick_label_complex(self, tick_label_theta, tick_label_radius_dict={'R': 1.15}):
        return maths.to_complex(tick_label_theta, tick_label_radius_dict['R'])


    def ideogram_complex(self, ideogram_coord_config, SUM, degreerange=[0,360], ideogram_radius_dict={'R0': 1, "R1": 1.1}):
        
        ## this will also be used for ring background!!!
        

        ideogram_theta_list = self.ideogram_theta_list(ideogram_coord_config, SUM, degreerange)
        radius_list = [*map(lambda x: np.concatenate((np.ones(len(x)//2)*ideogram_radius_dict['R0'],
                                                      np.ones(len(x)//2)*ideogram_radius_dict['R1'])), ideogram_theta_list
                            )]
        ideogram_complex_list = [*map(lambda x, y: maths.to_complex(x, y), ideogram_theta_list, radius_list)]

        #ideogram_complex_list_0 = maths.to_complex(self.ideogram_theta_list(SUM, degreerange), ideogram_radius_dict['R0'])
        #ideogram_complex_list_1 = maths.to_complex([*map(lambda x: x[::-1], self.ideogram_theta_list(SUM, degreerange))], ideogram_radius_dict['R1'])
        #ideogram_complex_list = [*map(lambda x, y: np.concatenate((x, y)), ideogram_complex_list_0, ideogram_complex_list_1)]
        return ideogram_complex_list


    def ideogram_path(self, ideogram_complex):
        # works for both ring and ideogram!
        ideogram_path_array_list = [*map(lambda x: np.column_stack((np.concatenate((np.full(1, 'M'), np.full(len(x)-1, 'L'))), x.real, x.imag)).ravel(), ideogram_complex)]
        ideogram_path_string_list = [*map(lambda x: " ".join(x) + ' Z', ideogram_path_array_list)]
        return ideogram_path_string_list

    def tick_path(self, tick_complex):
        path_array_list = [*map(lambda x: np.column_stack((np.tile(['M', 'L'], len(x)//2), x.real, x.imag)).ravel(), tick_complex)]
        path_string_list = [*map(lambda x: " ".join(x), path_array_list)]
        path_string = " ".join(path_string_list)
        return path_string


    def data_complex(self,
                     ideogram_coord_config, 
                     data_array,
                     category,
                     radius_dict,
                     SUM,
                     custom_offset_degree=False,
                     degreerange=[0,360],
                     return_path=True
                     ):

        # return_path=False: Traces for area plot, we don't add the return arc which would have ruined the hovertext arrangement
        # return_path=True: Layout for area plot, we need the return arc to close the area so that fill color can be applied

        if category == 'annotation':
            constant = np.pi/180
            data_theta = maths.to_theta(data_array[:,1], SUM, degreerange=degreerange)
            if custom_offset_degree:
                data_theta += data_array[:,3]*constant
            Data_complex = maths.to_complex(data_theta, radius_dict['R'])
            return Data_complex

        elif category in ['histogram', 'heatmap', 'cytoband']:
            # chr_name start end val



            assert self.ideogram_theta_list(ideogram_coord_config, SUM, degreerange=degreerange)
            data_theta_interval = maths.to_theta(data_array[:,1:3], SUM, degreerange=degreerange)
            data_complex_0 = maths.to_arc(data_theta_interval, radius_dict['R0'], self.ideogram_theta_list(ideogram_coord_config, SUM, degreerange=degreerange))

            if category == 'histogram':
                # when radius is dependent upon data_val
                data_radius = maths.val2radius(data_array[:,3], radius_dict['R0'], radius_dict['R1'])

                data_complex_1 = [*map(lambda x: x[::-1], maths.to_arc(data_theta_interval, data_radius, self.ideogram_theta_list(ideogram_coord_config, SUM, degreerange=degreerange)))]

            else:
                data_complex_1 = [*map(lambda x: x[::-1], maths.to_arc(data_theta_interval, radius_dict['R1'], self.ideogram_theta_list(ideogram_coord_config, SUM, degreerange=degreerange)))]
 

            Data_complex = [*map(lambda x, y: np.concatenate((x,y)), data_complex_0, data_complex_1)]
            return Data_complex

        elif category == 'scatter':
            data_radius = maths.val2radius(data_array[:,2], radius_dict['R0'], radius_dict['R1'])
            data_theta = maths.to_theta(data_array[:,1], SUM, degreerange=degreerange)
            Data_complex = maths.to_complex(data_theta, data_radius)
            return Data_complex

        elif category in ['line', 'area']:
            # chr_name pos val
            # get the unique chromosome index
            data_chr = data_array[:,0]
            _, indices = np.unique(data_chr, return_index=True)
            indices = np.sort(indices)

            data_radius = maths.val2radius(data_array[:,2], radius_dict['R0'], radius_dict['R1'])
            data_theta = maths.to_theta(data_array[:,1], SUM, degreerange=degreerange)
            data_complex_array = maths.to_complex(data_theta, data_radius)
            # splitting complex_array to a list by chromosome:
            Data_complex = np.split(data_complex_array, indices[1:])

            if category == 'area' and return_path is True:
                data_theta_list = np.split(data_theta, indices[1:])
                # theta interval should be the first and last theta value of each chromosome
                data_theta_interval = np.column_stack((np.array([*map(lambda x: x[0], data_theta_list)]), np.array([*map(lambda x: x[-1], data_theta_list)])))

                # the return arc_complex, needs to be reversed
                data_arc_complex = [*map(lambda x: x[::-1], maths.to_arc(data_theta_interval, radius_dict['R0'], self.ideogram_theta_list(ideogram_coord_config, SUM, degreerange)))]
                Data_complex = [*map(lambda x, y: np.concatenate((x,y)), Data_complex, data_arc_complex)]
            return Data_complex


        elif category == 'highlight':
            ### R0, R1 is defined from the file
            # chr_name start end R0 R1 opacity
            interval_theta_array = maths.to_theta(data_array[:,1:3], SUM, degreerange=degreerange)

            arc_complex_list_0 = maths.to_arc(interval_theta_array, radius_dict['R0'], self.ideogram_theta_list(ideogram_coord_config, SUM, degreerange=degreerange))


            arc_complex_list_1 = [*map(lambda x: x[::-1, ], maths.to_arc(interval_theta_array, radius_dict['R1'], self.ideogram_theta_list(ideogram_coord_config, SUM, degreerange)))]
            Data_complex = [*map(lambda x, y: np.concatenate((x,y)), arc_complex_list_0, arc_complex_list_1)]
            return Data_complex

        elif category == 'tile':
            # chr_name start end data_val
            interval_theta_array = maths.to_theta(data_array[:,1:3], SUM, degreerange=degreerange)
            data_radius = maths.val2radius(data_array[:,3], radius_dict['R0'], radius_dict['R1'])
            Data_complex = maths.to_arc(interval_theta_array, data_radius, self.ideogram_theta_list(ideogram_coord_config, SUM, degreerange=degreerange))
            return Data_complex

        elif category == 'connector':
            # chr_name pos0 pos1 

            ratio_constant = radius_dict['ratio'][2]/radius_dict['ratio'][0]
            data_theta = maths.to_theta(data_array[:,1:3], SUM, degreerange=degreerange)

            a0 = maths.to_complex(data_theta[:,0], radius_dict['R0'])
            a1 = maths.to_complex(data_theta[:,0], radius_dict['R0'] + (radius_dict['R1'] - radius_dict['R0'])*radius_dict['ratio'][0])
            a3 = maths.to_complex(data_theta[:,1], radius_dict['R1'])
            a2 = a3 + (a0-a1)*ratio_constant
            Data_complex = np.column_stack((a0, a1, a2, a3))
            return Data_complex

        elif category in ['link', 'ribbon', 'twistedribbon']:
            # R0 is bezier radius, R1 is ends radius
            # chr0 start0 end0 chr1 start1 end1
            assert radius_dict['R0'] is not None
            assert radius_dict['R1'] >= radius_dict['R0']
            data_coord_array = np.column_stack((data_array[:,1:3], data_array[:,4:6]))
            Data_complex = maths.bezier_complex(data_coord_array, radius_dict['R1'], radius_dict['R0'], SUM, type=category, degreerange=degreerange)
            return Data_complex

    def data_path(self,
                  ideogram_coord_config, 
                  category,
                  data_complex,
                  SUM,
                  degreerange=[0,360],
                  radius_dict={},
                  interval_theta_array_0=None,
                  interval_theta_array_1=None
                  ):
        length = len(data_complex)
    
        if category in ['histogram', 'heatmap', 'cytoband', 'tile', 'highlight', 'area', 'ring']:
            # data_complex is a ndarray list
            path_array_list = [*map(lambda x: np.column_stack((np.concatenate((np.full(1, 'M'), np.full(len(x)-1, 'L'))), x.real, x.imag)).ravel(), data_complex)]
            
            if category in ['line', 'tile']:
                # these paths dont form a loop shape, 
                path_string_list = [*map(lambda x: " ".join(x), path_array_list)]
            
            else:
                path_string_list = [*map(lambda x: " ".join(x) + ' Z', path_array_list)]
            return path_string_list
        
        elif category == 'connector':
            ## data_complex is a ndarray
            # in case where connector has customcolor, we need the path_string_list, otherwise we can just join
            path_array = np.column_stack((np.full(length, 'M'), data_complex[:,0].real, data_complex[:,0].imag, 
                                          np.full(length, 'L'), data_complex[:,1].real, data_complex[:,1].imag, 
                                          np.full(length, 'L'), data_complex[:,2].real, data_complex[:,2].imag, 
                                          np.full(length, 'L'), data_complex[:,3].real, data_complex[:,3].imag
                                          ))
            path_string_list = [*map(lambda x: " ".join(x), path_array)]
            return path_string_list


        elif category == 'link':

            if data_complex.ndim == 1:
                data_complex = data_complex.reshape((len(data_complex)//6, 6))
            length = len(data_complex)

            path_array = np.column_stack((np.full(length, 'M'), data_complex[:,0].real, data_complex[:,0].imag,
                                          np.full(length, 'Q'), data_complex[:,2].real, data_complex[:,2].imag,
                                          np.full(length, ''), data_complex[:,4].real, data_complex[:,4].imag,
                                          np.full(length, 'M'), data_complex[:,1].real, data_complex[:,1].imag,
                                          np.full(length, 'Q'), data_complex[:,3].real, data_complex[:,3].imag,
                                          np.full(length, ''), data_complex[:,5].real, data_complex[:,5].imag
                                          ))
            path_string_list = [*map(lambda x: " ".join(x), path_array)]


            # keep output pathstring as list in case we need custom color
            return path_string_list

        elif category in ['ribbon', 'twistedribbon']:
            ## Need to append two arcs, this is kinda complicated and I couldn't think of an easier way to do this
            if data_complex.ndim == 1:
                data_complex = data_complex.reshape((len(data_complex)//6, 6))
            length = len(data_complex)

            assert radius_dict['R1'] is not None
            assert interval_theta_array_0 is not None
            assert interval_theta_array_1 is not None

            arc0_complex_list = maths.to_arc(interval_theta_array_0, radius_dict['R1'], self.ideogram_theta_list(ideogram_coord_config, SUM, degreerange=degreerange))
            arc1_complex_list = maths.to_arc(interval_theta_array_1, radius_dict['R1'], self.ideogram_theta_list(ideogram_coord_config, SUM, degreerange=degreerange))

            if category == 'ribbon':
                link0_path_array = np.column_stack((np.full(length, 'M'), data_complex[:,5].real, data_complex[:,5].imag,
                                                    np.full(length, 'Q'), data_complex[:,2].real, data_complex[:,2].imag,
                                                    np.full(length, ''), data_complex[:,0].real, data_complex[:,0].imag
                                                    ))

                link1_path_array = np.column_stack((np.full(length, 'M'), data_complex[:,1].real, data_complex[:,1].imag,
                                                    np.full(length, 'Q'), data_complex[:,3].real, data_complex[:,3].imag,
                                                    np.full(length, ''), data_complex[:,4].real, data_complex[:,4].imag
                                                    ))
            else:
                ## reverse arc1 for twistedribbon!
                arc1_complex_list = [*map(lambda x: x[::-1], arc1_complex_list)]

                link0_path_array = np.column_stack((np.full(length, 'M'), data_complex[:,4].real, data_complex[:,4].imag,
                                                    np.full(length, 'Q'), data_complex[:,2].real, data_complex[:,2].imag,
                                                    np.full(length, ''), data_complex[:,0].real, data_complex[:,0].imag
                                                    ))

                link1_path_array = np.column_stack((np.full(length, 'M'), data_complex[:,1].real,data_complex[:,1].imag,
                                                    np.full(length, 'Q'), data_complex[:,3].real, data_complex[:,3].imag,
                                                    np.full(length, ''), data_complex[:,5].real, data_complex[:,5].imag
                                                    ))
                
            # because arc are always in the middle of the PATH, so only use L
            arc0_path_array_list = [*map(lambda x: np.column_stack((np.full(len(x), 'L'), x.real, x.imag)).ravel(), arc0_complex_list)]
            arc1_path_array_list = [*map(lambda x: np.column_stack((np.full(len(x), 'L'), x.real, x.imag)).ravel(), arc1_complex_list)]
            
            arc0_path_string_list = [*map(lambda x: " ".join(x), arc0_path_array_list)]
            arc1_path_string_list = [*map(lambda x: " ".join(x), arc1_path_array_list)]

            link0_path_string_list = [*map(lambda x: " ".join(x), link0_path_array)]
            link1_path_string_list = [*map(lambda x: " ".join(x), link1_path_array)]

            path_string_list = [*map(lambda a, b, c, d: " ".join([a, b, c, d]), link0_path_string_list, arc0_path_string_list, link1_path_string_list, arc1_path_string_list)]
            return path_string_list
        else:
            raise ValueError('Please select a supported category for path, scatter, line and annotation does not need path')
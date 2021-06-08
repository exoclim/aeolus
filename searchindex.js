Search.setIndex({docnames:["api/calc","api/const","api/coord","api/core","api/index","api/io","api/model","api/plot","api/region","api/subset","api/synthobs","changelog","contributing","examples/00_Constants","examples/01_Loading_Data","examples/02_Model_Field_Names","examples/03_Transmission_Spectrum","examples/index","index","install"],envversion:{"sphinx.domains.c":2,"sphinx.domains.changeset":1,"sphinx.domains.citation":1,"sphinx.domains.cpp":3,"sphinx.domains.index":1,"sphinx.domains.javascript":2,"sphinx.domains.math":2,"sphinx.domains.python":3,"sphinx.domains.rst":2,"sphinx.domains.std":2,"sphinx.ext.intersphinx":1,"sphinx.ext.viewcode":1,nbsphinx:3,sphinx:56},filenames:["api/calc.rst","api/const.rst","api/coord.rst","api/core.rst","api/index.rst","api/io.rst","api/model.rst","api/plot.rst","api/region.rst","api/subset.rst","api/synthobs.rst","changelog.rst","contributing.rst","examples/00_Constants.ipynb","examples/01_Loading_Data.ipynb","examples/02_Model_Field_Names.ipynb","examples/03_Transmission_Spectrum.ipynb","examples/index.rst","index.rst","install.rst"],objects:{"aeolus.calc":{calculus:[0,0,0,"-"],diag:[0,0,0,"-"],flux_h:[0,0,0,"-"],stats:[0,0,0,"-"],tl:[0,0,0,"-"]},"aeolus.calc.calculus":{d_dx:[0,1,1,""],d_dy:[0,1,1,""],d_dz:[0,1,1,""],deriv:[0,1,1,""],div_h:[0,1,1,""],integrate:[0,1,1,""]},"aeolus.calc.diag":{air_density:[0,1,1,""],air_potential_temperature:[0,1,1,""],air_temperature:[0,1,1,""],bond_albedo:[0,1,1,""],dry_lapse_rate:[0,1,1,""],flux:[0,1,1,""],geopotential_height:[0,1,1,""],ghe_norm:[0,1,1,""],heat_redist_eff:[0,1,1,""],horiz_wind_cmpnts:[0,1,1,""],meridional_mass_streamfunction:[0,1,1,""],precip_sum:[0,1,1,""],sfc_net_energy:[0,1,1,""],sfc_water_balance:[0,1,1,""],sigma_p:[0,1,1,""],toa_cloud_radiative_effect:[0,1,1,""],toa_eff_temp:[0,1,1,""],toa_net_energy:[0,1,1,""],water_path:[0,1,1,""],wind_speed:[0,1,1,""],zonal_mass_streamfunction:[0,1,1,""]},"aeolus.calc.flux_h":{horizontal_fluxes_through_region_boundaries:[0,1,1,""],net_horizontal_flux_to_region:[0,1,1,""]},"aeolus.calc.stats":{abs_coord_mean:[0,1,1,""],cumsum:[0,1,1,""],last_n_day_mean:[0,1,1,""],meridional_mean:[0,1,1,""],minmaxdiff:[0,1,1,""],normalize_cube:[0,1,1,""],region_mean_diff:[0,1,1,""],spatial:[0,1,1,""],spatial_mean:[0,1,1,""],spatial_quartiles:[0,1,1,""],time_mean:[0,1,1,""],vertical_mean:[0,1,1,""],zonal_mean:[0,1,1,""]},"aeolus.calc.tl":{regrid_to_rotated_pole_coordinates:[0,1,1,""],regrid_to_tidally_locked_coordinates:[0,1,1,""],rotate_winds_to_tidally_locked_coordinates:[0,1,1,""]},"aeolus.const":{add_planet_conf_to_cubes:[1,1,1,""],get_planet_radius:[1,1,1,""],init_const:[1,1,1,""]},"aeolus.coord":{CoordContainer:[2,2,1,""],add_binned_coord:[2,1,1,""],add_cyclic_point_to_cube:[2,1,1,""],add_planet_calendar:[2,1,1,""],area_weights_cube:[2,1,1,""],check_coords:[2,1,1,""],coarsen_cube:[2,1,1,""],coord_delta_to_cube:[2,1,1,""],coord_to_cube:[2,1,1,""],ensure_bounds:[2,1,1,""],get_cube_datetimes:[2,1,1,""],get_cube_rel_days:[2,1,1,""],get_dim_coord:[2,1,1,""],get_xy_coords:[2,1,1,""],interp_cube_from_height_to_pressure_levels:[2,1,1,""],interp_cubelist_from_height_to_pressure_levels:[2,1,1,""],interp_to_cube_time:[2,1,1,""],isel:[2,1,1,""],nearest_coord_value:[2,1,1,""],not_equal_coord_axes:[2,1,1,""],regrid_3d:[2,1,1,""],replace_z_coord:[2,1,1,""],roll_cube_0_360:[2,1,1,""],roll_cube_pm180:[2,1,1,""],vertical_cross_section_area:[2,1,1,""],volume_weights_cube:[2,1,1,""]},"aeolus.coord.CoordContainer":{__init__:[2,3,1,""]},"aeolus.core":{AtmoSim:[3,2,1,""],AtmoSimBase:[3,2,1,""],Run:[3,2,1,""]},"aeolus.core.AtmoSim":{sfc_water_balance:[3,4,1,""],sigma_p:[3,4,1,""],toa_net_energy:[3,4,1,""],wind_speed:[3,4,1,""]},"aeolus.core.AtmoSimBase":{"const":[3,4,1,""],__init__:[3,3,1,""],description:[3,4,1,""],extract:[3,3,1,""],from_parent_class:[3,3,1,""],model:[3,4,1,""],name:[3,4,1,""]},"aeolus.core.Run":{"const":[3,4,1,""],__init__:[3,3,1,""],add_data:[3,3,1,""],load_data:[3,3,1,""],model:[3,4,1,""],name:[3,4,1,""],proc_data:[3,3,1,""],to_file:[3,3,1,""]},"aeolus.io":{load_data:[5,1,1,""],load_multidir:[5,1,1,""],load_vert_lev:[5,1,1,""],save_cubelist:[5,1,1,""]},"aeolus.model":{um:[6,4,1,""],um_stash:[6,4,1,""]},"aeolus.plot":{cart:[7,0,0,"-"],mpl:[7,0,0,"-"],pv:[7,0,0,"-"],text:[7,0,0,"-"]},"aeolus.plot.cart":{GeoAxesGrid:[7,2,1,""],label_global_map_gridlines:[7,1,1,""]},"aeolus.plot.cart.GeoAxesGrid":{__init__:[7,3,1,""]},"aeolus.plot.mpl":{MidpointNormalize:[7,2,1,""],add_custom_legend:[7,1,1,""]},"aeolus.plot.mpl.MidpointNormalize":{__init__:[7,3,1,""]},"aeolus.plot.pv":{grid_for_scalar_cube_sph:[7,1,1,""],grid_for_vector_cubes_sph:[7,1,1,""]},"aeolus.plot.text":{fmt_lonlat:[7,1,1,""],subplot_label_generator:[7,1,1,""]},"aeolus.region":{Region:[8,2,1,""]},"aeolus.region.Region":{__init__:[8,3,1,""],add_to_ax:[8,3,1,""],constraint:[8,4,1,""],description:[8,4,1,""],from_cube:[8,3,1,""],name:[8,4,1,""]},"aeolus.subset":{DimConstr:[9,2,1,""],extract_last_month:[9,1,1,""],extract_last_n_days:[9,1,1,""],l_range_constr:[9,1,1,""],unique_cubes:[9,1,1,""]},"aeolus.subset.DimConstr":{__init__:[9,3,1,""]},"aeolus.synthobs":{calc_stellar_flux:[10,1,1,""],calc_transmission_spectrum_day_night_average:[10,1,1,""],read_normalized_stellar_flux:[10,1,1,""],read_spectral_bands:[10,1,1,""]},aeolus:{"const":[1,0,0,"-"],coord:[2,0,0,"-"],core:[3,0,0,"-"],io:[5,0,0,"-"],region:[8,0,0,"-"],subset:[9,0,0,"-"],synthobs:[10,0,0,"-"]}},objnames:{"0":["py","module","Python module"],"1":["py","function","Python function"],"2":["py","class","Python class"],"3":["py","method","Python method"],"4":["py","attribute","Python attribute"]},objtypes:{"0":"py:module","1":"py:function","2":"py:class","3":"py:method","4":"py:attribute"},terms:{"0":[0,1,2,7,13,14,16,18,19],"00":14,"00025":16,"00033333":16,"0005":16,"001":16,"01":16,"02":11,"03":11,"04":11,"05":[11,14],"05797807":13,"05887":16,"06":14,"07":[11,14],"08":11,"09":14,"1":[0,2,7,10,13,14,15,16],"10":[0,7,13,14,16],"100":16,"10000":16,"10001":16,"11":[13,14,16],"12":[11,13,14,16],"120":16,"121":7,"123":13,"1234":12,"13":[13,14],"14":13,"1400":16,"144":14,"15":[13,16],"16":2,"17":11,"180":[2,8],"1808":16,"194":16,"1994":0,"1d":[0,2],"1e":14,"1e6":16,"2":[0,1,2,3,7,10,13,14,15,16],"20":16,"2004":14,"2008":14,"2013":0,"2015":0,"2017":0,"2018":16,"2019":11,"2020":[11,18],"2021":11,"21":[0,11,14],"25":[7,11],"2501000":13,"25w":7,"27195379":16,"28":11,"287":13,"2d":7,"3":[13,14,15,16,18,19],"30":[0,11],"300":16,"31":11,"35":16,"360":2,"36000":14,"36360":14,"365":[0,9],"36720":14,"39":[13,14,15,16],"3d":[7,18],"4":[0,13,14,15,16],"40":16,"475026500":16,"481":16,"5":[0,13,14,15,16],"500":16,"5078378":16,"6":[0,13,14,15,16],"6371229":1,"7":[13,14,15,16,19],"8":[13,14,15,16,19],"80":[0,2,3,7,8,9,10],"85538168664425":16,"89":7,"8e6":16,"9":[13,14,15,16,19],"90":14,"break":4,"case":[2,10,11,13],"class":[0,2,3,7,8,9,10,11,15],"const":[0,1,2,3,11,13,14],"default":[0,1,2,7,16],"do":[11,12,17,19],"float":[0,2,7,10],"function":[0,3,4,5,10,11,14,15,18],"import":[0,7,12,13,14,15,16,19],"int":[2,3,7],"long":14,"new":[0,2,3,4,7,11,12,19],"public":4,"return":[0,1,2,3,5,7,8,10,11,12],"short":[0,13,14],"switch":7,"true":[0,2,3,7],"try":4,A:[0,2,3,7,8,13],By:[0,2,11],For:[0,13,14,15],If:[0,1,2,3,4,7,8,12,19],In:[0,13],It:[4,12,13,14,15,18],The:[0,2,3,7,8,10,12,13,14,15,16,18,19],Then:13,There:15,These:17,To:[12,14,16,19],__call__:7,__init__:[2,3,7,8,9],__path__:19,__repr__:11,_build:12,_cube:14,_subplot:7,a_:0,ab:[0,16],about:16,abov:[12,13],abs_coord_mean:[0,11],absolut:0,absorpt:16,academ:16,access:14,account:[12,16],action:11,activ:[2,4],ad:[2,13],adapt:[2,11],add:[0,1,2,3,7,8,11,12,14],add_binned_coord:[2,11],add_binned_lon_lat:11,add_custom_legend:[7,11],add_cyclic_point:2,add_cyclic_point_to_cub:2,add_data:3,add_planet_calendar:2,add_planet_conf_to_cub:[1,11],add_to_ax:8,addit:[3,5],aeolu:[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19],aeolus_data:[11,12,13,16],after:19,again:12,aggr:0,aggreg:0,air:[0,13],air_dens:[0,11],air_potential_temperatur:[0,11],air_temperatur:[0,11],al:[0,16,18],albedo:0,all:[0,2,11,12,13,19],allow:11,along:[0,2,7,11,13,14],alphabet:7,alreadi:19,also:[12,13,14,15,19],altern:12,altitud:[0,15],alwai:7,among:9,amount:15,an:[0,2,3,7,10,11,12,13,14,16],anaconda:19,analysi:[0,3,18],ani:[12,17,19],anoth:[0,11],antistellar:0,anybodi:19,anyth:3,api:[11,18],appear:12,append:[3,11],appli:[0,2,11],appropri:[0,3,13],april:11,ar:[0,1,2,3,4,7,11,12,13,14,15,16,17,19],arang:16,arbitrari:2,area:2,area_weights_cub:2,arg:[0,3],argument:[3,7,8],argumenterror:2,around:[7,11],arrai:[2,5,7,10,13,16],articl:16,arxiv:16,assign:3,associ:8,assum:[0,3],atmosflow:[2,11],atmosim:[3,11,14],atmosimbas:[3,11],atmospher:[0,2,3,10,13,14,16,18],atom:16,attach:11,attempt:[0,2,3],attribut:[0,1,2,3,5,11,14],au:[10,13,14],auto:2,autodoc:11,automat:[7,8,12],autoscale_non:7,aux_attr:5,auxiliari:[2,11],avail:[0,13,16],averag:[0,2,3,10,11,16],avogadro:[13,14],awar:13,ax:[2,7,8,16],axes_class:7,axes_grid1:7,axes_grid:7,axesgrid:7,axesgrid_kw:7,axessubplot:7,axi:[0,2,7],axis_weight:0,b:12,balanc:11,band:[10,16],bar:7,base:[2,3,7,8,9,11,15],basic:[11,13,14],becom:0,been:14,below:[12,14,16,19],best:17,better:[4,7,11,12],between:0,bin:2,black:12,blah:7,block:2,bold:0,boltzmann:[13,14],bond:0,bond_albedo:0,bool:[0,2,3,7,8],both:19,bottom:7,bound:[2,14],boundari:[7,8],branch:12,broadcast:[2,11],browser:12,bug:[11,12],build:[7,11,12,19],built:[7,19],c0:7,c1:7,c:[1,12,13,19],c_p:13,cach:[3,11],calc:[0,11],calc_stellar_flux:[10,16],calc_transmission_spectrum_day_night_averag:[10,16],calcul:[3,4,10,11,16],calculu:11,call:[7,14,15],callabl:[2,3],can:[8,12,13,14,15,19],cart:7,cartopi:[2,7,19],cd:[12,19],cdot:0,cell:14,cell_arrai:7,center:10,cfg:12,chang:[0,4,11,12],changelog:18,charact:3,check:[2,10,12,16],check_coord:[2,11],checkout:12,children:3,circul:18,classmethod:[3,8],clean:11,clear:0,climat:18,climatolog:0,clip:7,clone:12,close:12,closer:7,closest:2,cloud:0,cloud_wat:0,cluster:19,co:[0,2],coarsen_cub:2,code:[7,14,18],codebas:4,collaps:0,collect:1,color:7,colormap:7,column:16,com:[2,12,16,19],come:19,command:[12,19],commit:12,common:11,commonli:[0,15],compon:[0,3,7],comput:[2,3,12,19],concaten:16,conda:[11,12],condens:0,condensible_dens:[0,3,13,14],condensible_heat_vapor:[13,14],condensible_molecular_weight:13,condit:16,config:12,configur:[3,14],consist:2,const_dir:3,const_from_attr:11,constain:1,constant:[0,2,3,4,10,11,14,17],constcontain:[0,1,2,3],constraint:[0,3,4,8,11],contain:[0,1,2,3,5,9,11,14,15],contant:13,contourf:16,contributor:[11,18],conv:0,convective_rainfall_flux:14,convective_snowfall_flux:14,conveni:[11,12,13,14],convert:[0,2,7],coord:[0,2,9,11,16],coord_delta_to_cub:2,coord_nam:2,coord_to_cub:[2,11],coordcontain:2,coordin:[1,3,4,7,8,9,11,14,15],coordinatecollapseerror:11,copi:[2,12],core:[4,11,14],correct:[12,16],correctli:16,correspond:2,coupl:13,cov:12,coveragerc:12,cre:0,cre_:0,creat:[1,2,3,7,8,11,12,13,15],cross:2,ctrl:12,cube1:2,cube2:2,cube2d:2,cube:[0,1,3,4,5,7,8,10,11,13,14,16],cube_in:2,cube_src:2,cube_tgt:2,cubelist:[0,1,2,3,5,9,11,14],cumsum:0,cumul:0,current:0,custom:[7,13,15],cut:14,cv_rain:14,cwd:[13,14,16],cyclic:2,d:0,d_dx:0,d_dy:0,d_dz:0,dai:[0,1,2,9,13,14],data:[0,2,3,5,7,11,12,13,14,16],dataclass:[1,13],dataset:15,date:11,datetim:2,days_in_dai:2,days_in_month:2,days_in_year:2,daysid:[0,10,16],dc:9,deal:[2,3,5],decemb:11,decor:11,defeat:7,defin:[0,7,14,15],deg:13,degre:[2,7,8,13,14],delta:[2,11],densiti:0,depend:[11,19],deprec:11,deriv:[0,3,11,13],describ:10,descript:[3,8,14],desir:2,develop:[4,19],diag:[0,11],diagnost:[3,11],dict:[3,5,7],dictionari:[4,5,7],differ:[0,15],differenti:0,dim_constr_:11,dimconstr:[9,11],dimcoord:2,dimens:[0,2,7,9,11,14],dimension:[0,2,9],directli:[10,19],directori:[1,5,12,13,19],disk:5,distribut:19,div_h:0,diverg:[0,11],djf:2,doc:[2,12],docstr:11,document:[4,11],doe:[2,3,12,15],doesn:19,domain:[0,3,10,16],don:19,dot:12,download:[12,13],downward:[0,3],dp:0,drawn:7,drive:3,dry:[0,13],dry_air_gas_const:[0,13],dry_air_molecular_weight:[13,14],dry_air_spec_heat_press:[13,14],dry_lapse_r:[0,11],dt:0,dtype:16,dummi:13,dummyconst:13,duplic:[9,11,15],dynam:3,dz:0,e:[0,2,3,7,11,12],each:[11,13,14],earth:[0,1,2,3,11],earth_const:[0,13],earth_dai:[13,14],earthconst:[1,13],easier:3,east:7,east_bound:8,eastern:8,eccentr:[13,14],eff:0,effect:[0,10],effici:0,either:[14,19],element:2,els:19,emul:11,end:12,energi:[0,3,11],ensure_bound:2,env:12,environ:12,eq:0,equal:2,equat:11,equival:14,error:[11,12,19],et:[0,16,18],eta:0,evapor:[0,3],ex_const:13,exampl:[0,1,4,7,9,11,13,14,16,18],exclud:2,exist:[7,12],exner:0,exoclim:[12,13,16,19],exoclimatolog:16,expect:12,experi:14,extens:[1,2],extent:[10,16],extra:[14,19],extract:[0,3,9,11],extract_cub:14,extract_last_month:9,extract_last_n_dai:9,extrema:0,f4:16,f5:12,f:12,f_:[0,10],face:19,factor:7,failur:19,fall:7,fals:[0,2,3,7,8],familiar:12,feel:4,few:[11,13,17],field:[0,2,3,7,8,9,10,14,16],fig:[7,16],figsiz:16,figur:[7,18],file:[1,3,5,10,12,13,14,16],filenam:14,filter:11,filterwarn:16,find:[2,19],first:[0,3,7,11,13,16],fix:[11,12],flag:12,flake8:12,flux:[3,10,11,16],flux_h:0,fmt_lonlat:7,folder:[1,3,11,12],follow:[12,13,15,16,19],foo:7,forecast:2,forecast_period:14,forecast_reference_tim:14,forg:[11,19],forget:12,fork:12,format:[7,12],found:2,four:16,frac:[0,10],fraction:2,free:19,from:[0,1,2,3,5,7,8,9,10,11,12,13,14,15,16],from_cub:8,from_parent_class:3,full:5,func:3,func_arg:3,further:12,futur:4,g:[0,2,7,12,15],ga:13,gamma:0,gear:18,gen:13,gener:[0,1,2,7,11,18],generalis:2,geoax:7,geoaxesgrid:7,geoaxessubplot:7,geograph:[0,8],geopotenti:0,geopotential_height:[0,11],get:[0,1,2,3,4,19],get_cube_datetim:[2,11],get_cube_rel_dai:[2,11],get_dim_coord:2,get_planet_radiu:1,get_xy_coord:[2,11],ghe:0,ghe_norm:0,git:[12,19],github:[2,11,12,16],given:[0,1,2,3,7,8],global:[0,2,3,7],go:12,gone:19,graviti:[0,1,13,14],greenhous:0,grid:[0,2,3,7,11],grid_for_scalar_cube_sph:7,grid_for_vector_cubes_sph:7,gridlin:7,guess_coord_axi:2,guid:[11,18],h_max:9,h_min:9,ha:[12,13,15,18,19],had:19,hand:15,haqq:0,hartmann:0,have:[0,2,3,9,12,13,14,15,19],heat:[0,13],heat_redist_eff:0,heavili:18,height:[0,2,5,7,10,16],height_domain:16,helper:11,hemispher:7,here:[4,13,16,17],heurist:2,hierarchi:4,high_type_cloud_area_fract:14,highli:19,horiz_wind_cmpnt:[0,11],horizont:[2,11],horizontal_fluxes_through_region_boundari:0,hot:[10,16],hour:14,how:[16,17],html:12,http:[2,12,16,19],i:[0,3,4,7],i_cub:0,ic:0,ice_wat:0,ident:1,identif:3,identifi:7,idx:2,ignor:16,illustr:17,imagegrid:7,improv:[4,11,12],inact:2,includ:[2,7,12],include_pressur:2,independ:13,index:[2,16,18,19],indic:15,individu:14,info:14,inform:[3,16,19],inherit:3,init_const:[0,1,3,13],initi:7,initialis:[3,7,9,11],inp_data:14,input:[0,2,3,5,7,11,16],input_cub:0,instal:[12,18],instanc:[0,7],instanti:[2,3,8,14],instead:[7,12,15],institut:19,instruct:16,int_:0,integ:7,integr:[0,11],interfac:11,intern:3,interp_all_to_pres_lev:11,interp_cube_from_height_to_pressure_level:[2,11],interp_cubelist_from_height_to_pressure_level:[2,11],interp_to_cube_tim:[2,11],interp_to_pres_lev:11,interp_to_single_pres_lev:11,interpol:[0,2,11],io:[5,11],iri:[0,1,2,3,5,7,8,9,10,11,13,14,16,18,19],isel:[2,11],issu:[12,19],item:10,iter:3,its:[0,1,3,9,11,19],j:[0,13,14],j_cube:0,januari:11,jja:2,json:[1,3,13,14],june:11,jupit:[10,16],jupyt:11,just:4,k:[13,14],keep:[14,16],kei:[7,14],keyword:[3,7,8,11],kg:[13,14],kind:0,km:9,know:19,kopparapu:0,kwarg:8,l:13,l_range_constr:9,label:[3,5,7,19],label_global_map_gridlin:7,label_nam:5,lam:3,lambda:0,languag:19,laps:0,last:[0,9],last_n_day_mean:0,lat:[7,10],lat_band_constr:0,lat_bin:2,latent:0,latest:19,latitud:[0,2,3,7,8,11,14,15],latlon23:19,lead:7,lecont:0,left:[0,7],leg_kw:7,legend:7,length:[2,9],letter:2,lev:5,lev_typ:5,level:[2,3,5,7,9,11,12,19],level_height:[3,7,9,15],leverag:18,lgpl:18,librari:[11,12,18],like:[2,3,4,7,12,13],limb:16,limit:8,line:[7,12,16],linearli:2,linux:19,liquid:0,liquid_wat:0,list:[0,1,2,3,5,10,14],lmd:15,lmdg:15,load:[3,5,11,13,14],load_cub:16,load_data:[3,5,11],load_multidir:5,load_vert_lev:5,loc:7,local:12,localhost:12,locat:10,lock:[11,14],log:16,login:16,lon:[7,10],lon_bin:2,lon_or_lat:7,long_nam:16,longitud:[0,2,3,7,8,14,15],longwav:0,look:[2,13],lost:0,low_type_cloud_area_fract:14,lower_wavelength_limit:[10,16],lsaffin:2,lt:16,lw:0,m01s00i002:15,m01s00i003:15,m01s00i150:15,m01s01i755:[10,16],m01s05i205:14,m2:[13,14],m:[0,1,2,10,12,13,14,16],mac:19,made:[0,2,3],magnitud:[0,3],mai:11,main:[3,19],make:[0,9,11,12],mam:2,manual:12,map:7,march:11,margin:8,margin_unit:8,maria:11,marker:7,mask:7,mass:0,master:12,match:2,mathtext:16,matplotlib:[8,11,16,19],matrix:11,max:0,maximum:[0,7],mean:[0,14],medium_type_cloud_area_fract:14,merg:[5,11],meridion:[0,7,11],meridional_mass_streamfunct:0,meridional_mean:0,messag:19,met:[10,14,15,16,18],meta:11,metadata:[0,2,3,11,14],method:[11,14,17],methodolog:16,metpi:[11,19],metr:[1,2],midpoint:7,midpointnorm:7,min:0,miniconda:19,minimum:[0,7],minmaxdiff:0,minor:[11,19],minu:[0,3],misc:11,misra:0,miss:12,mix:0,mm:0,mmsf:0,mnra:16,mode:[12,19],model:[0,2,3,4,7,8,9,10,11,13,16,17,18],model_typ:3,modul:[11,16,18],mol:[13,14],molar_gas_const:[13,14],molecul:16,month:[2,9],more:[2,11,12,14],move:11,mpl:7,mpl_toolkit:7,mu:16,multipl:5,multipli:0,must:[0,2,3],my:13,my_const:13,my_data_cub:0,my_dict:7,my_run:14,myx:9,mzsf:0,n:[0,2,7,9],nabla:0,name:[0,1,2,3,4,7,8,9,10,11,13,14,17],navig:12,nc:16,ncol:7,ndarrai:[2,10],nearest:2,nearest_coord_valu:2,necessari:[0,2,12],necessarili:17,need:[12,16],neg:7,net:0,net_horizontal_flux_to_region:0,netcdf:[11,16],new_branch_for_cool_addit:12,night:0,nightli:19,nightsid:[0,10,16],non:[0,2,11],none:[0,1,2,3,7,8,19],normal:[0,2,7,10,16,19],normalis:[0,7],normalize_cub:[0,11],normalized_stellar_flux:16,north:[0,7],north_bound:8,northern:8,not_equal_coord_ax:2,note:[0,2,7,10,12,13,14,16],notebook:11,notfounderror:2,novemb:11,now:11,np:[0,16],nrow:7,nrows_ncol:7,nu:10,number:[2,8,13],numpi:[2,5,10,16,19],o:[4,7],obj:[2,3],object:[2,3,5,7,8,9,13,14],obliqu:[13,14],observ:[4,17],obtain:16,occur:4,octob:11,off:[4,7],offic:[10,14,15,16,18],offset:7,often:12,olr:0,olr_:0,onc:19,one:[0,2,7,11,12,13,14,15],onli:[0,1,7,9,13,16],open:[12,19],oper:[0,13],option:[0,1,2,3,5,7,8,9,10,11],orbit:11,ordin:2,org:16,origin:[0,2],osr_:0,other:[0,9,11],otherwis:7,oup:16,out:11,output:[2,3,5,10,16,17,18],outsid:7,over:[0,7,8,10,11,16],overal:4,overlin:0,overrid:11,own:12,p:[0,3,10,11,15],p_:0,p_ref_frac:2,pa:[13,14],packag:[3,4,11,16,18,19],page:[11,12,18,19],paper:12,paramet:[0,1,2,3,5,7,8,9,10,11,16],parent:3,part:19,partial:0,particular:17,pass:[3,7,12,15],patch:8,path:[0,1,3,5,10,12,13,14,16],path_mask:5,path_to_fil:5,pathlib:[1,3,5,10,13,14,16],pep8:12,per:[0,10,16],period:2,permiss:19,phase:0,phi:0,php:16,physic:[0,3,4,14,17],pi:0,pip:[12,19],plane:0,planet:[0,1,2,3,11,13,14,16],planet_conf:[5,14],planet_domain_height:16,planet_radiu:16,planet_rotation_r:11,planet_top_of_atmospher:[10,16],planet_transmission_dai:[10,16],planet_transmission_night:[10,16],planetari:[1,10,16,18],pleas:[12,16,19],plot:[4,11,16],plt:[7,16],point:[0,2,7,8,16],point_arrai:7,pointer:3,pole:0,pole_lat:0,pole_lon:0,port:12,posit:7,possibl:[13,14],post:3,potenti:0,pp:14,practic:17,pre:[2,11,15],precip_sum:0,precipit:[0,3,11],prepar:[7,11],present:[0,16],pressur:[0,2,3,11],primari:19,primarili:18,print:[14,16,19],problem:19,proc:3,proc_data:3,procedur:7,process:[3,7,11,14],profil:12,project:[7,15],prompt:12,properti:[3,7,11],provid:[4,13,14,15,16,17],psi_m:0,psi_z:0,ptype:0,publish:16,pull:12,purpos:7,push:12,put:7,pv:[7,11],py:[12,19],pyplot:[7,16],pytest:[12,19],python:[11,12,18,19],pyvista:[11,19],q:0,qplt:16,quantiti:[0,3],quartil:0,quickplot:16,r:0,r_:[10,16],r_d:13,r_p:[10,16],r_planet:[0,2],radi:0,radiat:0,radiu:[0,1,2,10,11,13,14,16],rain:0,rais:[2,11],random:12,rang:[0,7,9],rate:0,ratio:[0,10],rcparam:16,re:19,read:[5,10,13,14,16],read_normalized_stellar_flux:[10,16],read_spectral_band:[10,16],real:0,recip:11,recommend:14,rect:7,rectangl:8,rectangular:[0,8],redistribut:0,reduc:0,refactor:11,refer:[0,2,11,12,18,19],reference_surface_pressur:[0,13,14],refresh:12,region:[0,2,4],region_a:0,region_b:0,region_mean_diff:0,regrid3d:11,regrid:[0,2,11],regrid_3d:2,regrid_to_rotated_pole_coordin:0,regrid_to_tidally_locked_coordin:0,regular:[3,16],regularli:2,rel:[2,12],relat:2,relax:9,releas:[11,18],relev:[0,1,2,3,7,8,9,10,12,14,16],relevel:2,reli:18,remain:7,remot:12,remov:[9,11],renam:[11,13],replac:[1,2,11],replace_z_coord:2,repo:19,report:[12,16],repositori:[11,13,18],repres:[2,7,17],represent:11,request:[2,12],requir:[2,16,19],respect:[0,2,7,8],rest:[12,13],restructur:11,result:[0,19],retriev:[0,2,3,11],review:12,rewrit:11,rho:[0,5],right:[0,12],roll:2,roll_cube_0_360:2,roll_cube_pm180:2,root:19,rotat:[0,11],rotate_winds_to_tidally_locked_coordin:0,rp_eff_over_r:16,rule:0,run:[3,5,11,12,19],run_start_dai:2,s:[0,1,2,3,4,7,8,11,13,14,16,18,19],s_:0,same:[0,2,7,15,16,19],sampl:14,sample_fil:14,sample_t1e_2d_mean:14,saturn:11,save:[3,5,11],save_cubelist:5,scalar3d:7,scalar:[0,1,3,8,13,14],scalar_cub:0,scalarcub:11,scale:7,scienc:4,scientif:17,search:[18,19],season:2,second:0,section:[2,4],see:[0,2,12,14],select:2,selector:2,self:3,semi_major_axi:[13,14],send:19,sensibl:0,separ:[11,12],septemb:11,sequenc:[2,7],sergeev:18,serv:13,server:12,session:12,set:[1,2,7],set_titl:16,set_xlabel:16,set_xscal:16,set_xtick:16,set_xticklabel:16,set_ylabel:16,setup:[12,19],sfc:0,sfc_net_energi:0,sfc_water_bal:[0,3,11],shape:[2,11],shift:8,shift_lon:8,shortcut:[0,11,15],shortwav:0,should:[1,7,12],show:12,shown:13,sigma:0,sigma_p:[0,3],silent:7,similar:3,simpl:11,simplic:11,simul:[2,3,14],singl:[2,3,14,16],skip:[2,12],skip_not_found:2,sky:0,slice:[2,9],small:12,smaller:8,snow:0,so:[7,12,13],socrat:[10,16],softwar:17,solar_const:[0,1,13,14],sole:3,some:[0,2,4,11,13,14,15],someth:19,son:2,sourc:[0,1,2,3,5,7,8,9,10,14,19],south:[0,7],south_bound:8,southern:8,sp:0,sp_sw_500ir_bd_hatp11:16,space:[0,2],span:2,spatial:[0,11],spatial_mean:[0,11],spatial_quartil:0,specif:[1,3,4,11,15],specifi:[0,2,14],spectra:16,spectral:[10,16],spectral_band_cent:16,spectral_band_index:[10,16],spectral_fil:[10,16],spectrum:[10,11],speed:[0,3],spheric:[0,3,7,11],sqrt:[0,3,10],src:11,standard:[7,11,15,19],start:[2,4,12],stash:[10,14,15],stat:0,state:10,stefan_boltzmann:[13,14],stellar:[10,16],stellar_constant_at_1_au:[10,16],stellar_flux:16,stellar_radiu:[10,16],step:[3,12],store:[3,11,13,14,15,16],str:[0,1,2,3,5,7,8,14,16],stra:0,stratifi:[2,11,19],stratiform_rainfall_flux:14,stratiform_snowfall_flux:14,streamfunct:[0,11],strict:9,stride:7,string:7,structur:10,structuredgrid:7,studi:17,style:12,styles_and_label:7,sub:[12,13],subclass:7,submodul:[3,11,13,15],subplot:[7,16],subplot_label_gener:7,subset:[3,4,8,11],substellar:0,suit:19,sum:[0,10,11,16],sum_:10,superrotation_index:11,support:19,sure:12,surfac:[0,2,3,11],surpris:7,sw:0,symbol:7,symmetr:11,synthet:[4,11,17],synthob:[10,16],system:1,t1e_exampl:14,t:[0,2,15,19],t_:0,tabl:15,take:[2,3,7,12,16],target:2,team:12,technic:11,temperatur:0,term:[12,18],test:[11,16,19],test_data:[13,14,16],tex:7,text:11,text_kw:7,textbook:12,than:8,thei:7,them:[0,2,3,7,12],therefor:7,theta:[0,5],thi:[0,3,7,12,13,14],those:11,through:[0,19],tick:7,tidal:[11,14],time:[0,2,3,9,11,14,15],time_coord:2,time_mean:[0,11],timestep:[3,11],titan:[2,11],titl:[7,16],tl:0,to_fil:3,toa:[0,3,10,11],toa_cloud_radiative_effect:0,toa_eff_temp:0,toa_net_energi:[0,3],top:[0,12,19],total:[0,3,10],toward:18,transform:0,transmiss:[10,11],transmission_spectra:16,transmit:10,trap1:14,trap1econst:14,trapezoid:0,trappist:14,travi:11,travisci:11,treat:2,tupl:[2,10],two:[0,11,16],type:[0,1,2,3,5,7,8,10,12,13,14,17],typo:[11,12],u4:16,u:[0,1,3,7,15],uk:[15,18],um:[0,5,6,11,14,15,16],um_:11,um_stash:[6,15],um_vers:14,under:[4,7,18],unequ:2,unifi:[10,14,16,18],unique_cub:[9,11],unit:[0,2,7,8,9,10,11,13,16],univers:13,up:[0,2,11],updat:[11,12],update_metadata:11,upon:18,upper_wavelength_limit:[10,16],upward_air_veloc:15,url:12,us:[0,2,3,4,7,8,12,13,14,15,16,17,18,19],usag:[14,17],user:3,usual:[0,10,16],util:[2,7,11],v3:11,v:[0,3,7,12,15],val:2,valli:0,valu:[0,2,7,13],vapor:13,vapour:0,variabl:[0,2,3,4,7,8,10,11,17],variou:[15,19],vec:0,vector3d:7,vector:[0,3,7],vector_scal:7,veloc:1,version:12,vert_coord:3,vertic:[0,2,3,5,7,11],vertical_constraint:0,vertical_cross_section_area:2,vertical_mean:[0,11],vertical_sum:11,vertical_wind_scal:7,via:[14,19],visibl:7,visualis:18,vmax:7,vmin:7,volum:[2,11],volume_weights_cub:2,w:[0,1,3,5,7,10,13,14,15,16],wa:[12,19],wai:[13,14,17],wait:12,want:19,warn:[12,16],warn_thresh:0,water:[0,3,13],water_heat_vapor:[13,14],water_molecular_weight:[13,14],water_path:0,water_vapour:0,wavelength:16,we:[12,13,19],weight:[0,2],weight_bi:0,weight_by_dens:0,welcom:12,west:7,west_bound:8,western:8,what:[5,12],when:2,where:[0,2,7,10,12,19],which:[2,10,12,13,14,15,16,18,19],whichev:7,whole:15,width:7,wiki:16,wildcard:3,wind:[0,3,7],wind_spe:[0,3,11],window:19,within:[2,7,9,14],without:1,work:[0,2,11,12,17,18,19],workflow:4,wp:0,wrap:5,wrapper:2,write:5,wrong:19,www:19,x:[0,2,7,9,14,15],x_i:0,x_wind:15,xarrai:[11,19],xlabel:16,xoff:7,xstride:7,xtick:[7,16],y:[0,2,7,9,15],y_wind:15,year:2,ylabel:16,yml:12,yoff:7,you:[4,12,14,16,19],your:12,your_user_nam:12,yourself:12,ystride:7,ytick:7,yx:9,z:[0,2,3,7,15],z_:0,z_i:0,z_offset:7,z_scale:7,zamyatina:11,zonal:[0,1,7,11],zonal_mass_streamfunct:0,zonal_mean:0},titles:["Science calculations","Physical constants","Cube coordinate functionality","Core hierarchy","API reference","I/O","Model-specific dictionaries of variable and coordinate names","Plotting functions","Regions","Subset cubes using constraints","Synthetic observations","Changelog","Contributor\u2019s Guide","Physical constants","Working with model output","Model variable names","Working with synthetic observations","Examples","Overview","Installation"],titleterms:{"0":11,"1":11,"10":11,"11":11,"2":11,"3":11,"4":11,"5":11,"6":11,"7":11,"8":11,"9":11,"function":[2,7],also:18,altern:19,api:4,calcul:0,calculu:0,cartographi:7,changelog:11,clone:19,code:12,conda:19,constant:[1,13],constraint:9,contain:13,content:18,contribut:12,contributor:12,coordin:[0,2,6],core:3,cube:[2,9],develop:12,diagnost:0,dictionari:6,document:12,earth:13,exampl:[15,17],extend:13,flux:0,from:19,gener:13,github:19,guid:12,helper:7,hierarchi:3,horizont:0,i:5,indic:18,instal:19,interfac:7,licens:18,lock:0,matplotlib:7,method:19,model:[6,14,15],name:[6,15],o:5,observ:[10,16],output:14,overview:18,physic:[1,13],plot:7,proper:19,pypi:19,pyvista:7,recommend:19,refer:4,region:8,relat:7,repositori:12,s:12,scienc:0,see:18,set:12,specif:6,spectrum:16,statist:0,subset:9,synthet:[10,16],tabl:18,test:12,text:7,tidal:0,transmiss:16,troubleshoot:19,unifi:15,up:12,us:9,variabl:[6,15],verifi:19,work:[14,16]}})
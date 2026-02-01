# This CASA pipescript is meant for use with CASA 6.6.1 and pipeline 2024.1.1.22
context = h_init()
context.set_state('ProjectSummary', 'observatory', 'Karl G. Jansky Very Large Array')
context.set_state('ProjectSummary', 'telescope', 'EVLA')
try:
    hifv_importdata(vis=['myCaldMS.ms'], datacolumns={'data': 'raw','corrected': 'regcal_contline_all'})
    hifv_flagtargetsdata()
    hifv_mstransform()
    hif_checkproductsize(maximsize=16384)
    hif_makeimlist(specmode='cont', datatype='regcal')
    hif_makeimages(hm_cyclefactor=3.0)
    hif_selfcal()
    hif_makeimlist(specmode='cont', datatype='selfcal')
    hif_makeimages(hm_cyclefactor=3.0)    
    hifv_pbcor()
    hifv_exportdata(imaging_products_only=True)
finally:
    h_save()

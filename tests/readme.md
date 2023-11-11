

# Parameters validation

## First validation and default parameter fitting

Directory: `1_validate_params`:

Files:
* First round of validation and default parameter set fitting
    * `utils_1.py` - file to be importes in other files (with param sets and K sets, names of files);
    * `1.0_iterate_params.py` - script to iterate over params and run simulation;
    * `1.1_iterate_results.py` - script to iterate over results of simulation, clump SNPs and summarize results. 
    * `1.2_check_parameters.ipynb` -- vizualization of these results.
* Second round of the revalidation of selected sets, same 4 files:
    * `utils_2.py` - file to be importes in other files (with param sets and K sets, names of files);
    * `2.0_iterate_best_params.py` - script to iterate over params and run simulation;
    * `2.1_iterate_best_results.py` - script to iterate over results of simulation, clump SNPs and summarize results. 
    * `2.2_validate_parameters.ipynb` -- vizualization of these results.
* ...


<!-- Other files which were used, but not included due to size: 
* `gencode.v37.annotation.gtf` -


Доисать: как получил значимые снипы (инпутфайлы для симуляции) когда симулировал ukb & finngen
И добавить к скрипту отрисовки актиуальные данные


Не забыть переименовать img/imgs -> image, in/out_data -> data



подгрузить данные и скрипт со второго сервака и картинки отрисовка и сами картинки с компа -->

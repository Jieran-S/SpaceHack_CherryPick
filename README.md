# SpaceHack - contributing modules (data, methods, metrics)

Our workflow is setup to allow everyone to contribute "modules" in their preference programming language (.. as long as that is either R or python). A module can either be a dataset, a computational method, or an evaluation metric.
![image](https://github.com/SpatialHackathon/SpaceHack2023/assets/114547/7c002916-0a90-4fe7-8745-489313bc0192)

This repository contains some templates and examples of how to implement your module so that they interface seamlessly with other modules in the workflow. For example, if you want to implement a new method, you do not need to worry about input data or evaluation metrics as long as you follow the template for reading input and writing output - if you correctly adhere to the input and output guidelines, you should be able to interface with our default data modules and default evaluation metrics modules. The default modules are:
 - data: LIBD Visium DLPFC dataset (4 samples, each with 3 replicates)
 - methods: BayesSpace and SpaGCN
 - evaluation metrics: ARI and V

### How to contribute a module

Module contribution will be managed via GitHub. The steps to contribute a module are:
 1. Create or claim a **GitHub issue** from the [SpaceHack issue board.](https://github.com/SpatialHackathon/SpaceHack2023/issues) that describes the module you want to implement. There are currently 90 issues to claim, but if you come up with a new idea, please **create** a new issue, add the appropriate **tags**, and **assign** the task to yourself.
 2. Add **metadata** to our metadata [spreadsheet](https://docs.google.com/spreadsheets/d/1QCeAF4yQG4bhZSGPQwwVBj_XF7ADY_2mK5xivAIfHsc/edit). Please fill in as much as you can as metadata is helpful! If you feel the need, please also add new columns or add additional notes. The metadata should be dded to the appropate tabs:
    - [datasets](https://docs.google.com/spreadsheets/d/1QCeAF4yQG4bhZSGPQwwVBj_XF7ADY_2mK5xivAIfHsc/edit#gid=1453488771)
    - [computational methods](https://docs.google.com/spreadsheets/d/1QCeAF4yQG4bhZSGPQwwVBj_XF7ADY_2mK5xivAIfHsc/edit#gid=0)
    - [evaluation metrics](https://docs.google.com/spreadsheets/d/1QCeAF4yQG4bhZSGPQwwVBj_XF7ADY_2mK5xivAIfHsc/edit#gid=4776337)
    - [simulations and technical evaluation](https://docs.google.com/spreadsheets/d/1QCeAF4yQG4bhZSGPQwwVBj_XF7ADY_2mK5xivAIfHsc/edit#gid=640974611)
 4. Now you are ready to create a new git **branch**
 5. Modify code
 6. Test
 7. Merge request
 8. Code review (?)

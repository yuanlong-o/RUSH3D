run("Trainable Weka Segmentation 3D");
wait(3000);
selectWindow("Trainable Weka Segmentation v3.2.33");
call("trainableSegmentation.Weka_Segmentation.loadClassifier", "D:\\RUSH3Dpipeline\\va.2\\segmentation_module\\classifier3d_20220303_1.model");
call("trainableSegmentation.Weka_Segmentation.getResult");
selectWindow("Classified image");
selectWindow("Trainable Weka Segmentation v3.2.33");
call("trainableSegmentation.Weka_Segmentation.getProbability");
selectWindow("Probability maps");
close("Trainable Weka Segmentation*");
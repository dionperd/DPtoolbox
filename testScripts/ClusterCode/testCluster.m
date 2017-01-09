function testCluster( filenum )


    load(['/home/mpib/perdikis/testClusterFolder/inputData/Input',filenum,'.mat'])

    out = 2*tool(Input);
    
    save(['/home/mpib/perdikis/testClusterFolder/Results/Res',filenum,'.mat'],'out')
    
end


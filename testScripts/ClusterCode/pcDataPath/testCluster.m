function testCluster( filenum )


    load(['/home/mpib/perdikis/testClusterFolder/Input',filenum,'.mat'])

    out = 2*tool(Input);
    
    save(['/home/mpib/perdikis/testClusterFolder/Output',filenum,'.mat'],'out')
    
end


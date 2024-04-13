for(const auto pop:population){
        std::cout<<"Poblacion:";
        for(size_t i=0; i<pop.size();i++){
            std::cout<<pop[i];
        }   
        std::cout<<std::endl;
    }

    for(auto allel:child1){
                std::cout<<allel;
            }
            std::cout<<"\n";

    std::vector<int> prueba={1,2,3,4,5,6,7,8,9};
    int pruebaV=calculateFitness(prueba);
    std::cout<<pruebaV<<std::endl;

child1=crossoverOX(population[selectedParents[i]], population[selectedParents[i+1]]);
                child1=mutationIn(child1);
                children.push_back(child1);
                child1.clear();
                child2=crossoverOX(population[selectedParents[i+1]], population[selectedParents[i]]);
                child2=mutationIn(child2);
                children.push_back(child2);
                child2.clear();
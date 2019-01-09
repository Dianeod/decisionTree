from pyclustering.utils import read_sample
from pyclustering.utils import timedcall

from pyclustering.samples.definitions import SIMPLE_SAMPLES
from pyclustering.samples.definitions import FCPS_SAMPLES

from pyclustering.cluster import cluster_visualizer
from pyclustering.cluster.cure import cure


def template_clustering(number_clusters, path, number_represent_points = 5, compression = 0.5, draw = True, ccore_flag = False):
    sample = read_sample(path)
    
    cure_instance = cure(sample, number_clusters, number_represent_points, compression, ccore_flag)
    (ticks, _) = timedcall(cure_instance.process)
    
    clusters = cure_instance.get_clusters()
    representors = cure_instance.get_representors()
    means = cure_instance.get_means()

    print("Sample: ", path, "\t\tExecution time: ", ticks, "\n")
    #print([len(cluster) for cluster in clusters])

    


template_clustering(2,'/home/mlteam/Documents/sample2.txt')



/*
 * Asessing the supertree clade by clade
 * input: Supertree and input trees in a newick file
 * output: a file containing usefull parameters
 * 
 * File:   main.cpp
 * Author: Ruchi Chaudhary
 *
 * Created on November 11, 2013, 2:04 PM
 */

#include <stdlib.h>
#include "argument.h"
#include "common.h"
#include "tree.h"
#include "tree_IO.h"
#include "tree_traversal.h"
#include "tree_LCA.h"
#include "tree_LCA_mapping.h"
#include "tree_name_map.h"
#include "tree_subtree_info.h"
#include <iomanip>

static const unsigned int NONODE = UINT_MAX;


/*
 * main file
 */
int main(int ac, char* av[]) {
    {
        std::ostringstream os; os << "command:";
        for (int i = 0; i < ac; i++) os << ' ' << av[i];
        MSG(os.str());
    }
    std::string trees_filename;
    std::string output_filename;
    unsigned int option = 0;
    bool taxonomy = false;
    {
        Argument a;
        a.add(ac, av);
        // help
        if (a.existArg2("-h","--help")) {
            MSG("options:");
            MSG("  -i [ --input ] arg      supertree and input trees (file in NEWICK format)");
            MSG("  -o [ --output ] arg     assessment");
            MSG(" -c 1 = stat, 2 = triplet count, 3 = RF distance, and 4 = Mark Wilkinson et al. based assess");
            MSG(" --tax");
            MSG("  -h [ --help ]           produce help message");
            MSG("");
            MSG("example:");
            MSG("  " << av[0] << " -i inputF.newick -o outputF.newick");
            exit(0);
        }

        // input trees
        if (a.existArgVal2("-i", "--input", trees_filename)) MSG("input file: " << trees_filename) else MSG("using standard input");
        
        // output file
        if (a.existArgVal2("-o", "--output", output_filename)) MSG("output file: " << output_filename);

        // argument option
        if (a.existArgVal("-c", option)) MSG("Option: " << option);
        if(option<1 || option >4)
            ERROR_exit("Not a valid option!");

        taxonomy = a.existArg("--tax");
        
        // unknown arguments?
        a.unusedArgsError();
    }

    //create output stream
    std::ofstream ouput_fs;
    if (!output_filename.empty()) {
        ouput_fs.open(output_filename.c_str());
        if (!ouput_fs) ERROR_exit("cannot write file '" << output_filename << "'");
    }
    std::ostream &output = output_filename.empty() ? std::cout : ouput_fs;

    // declaring variables
    aw::Tree s_tree;
    aw::idx2name s_taxa;
    std::vector<aw::Tree> g_trees;    
    std::vector<aw::idx2name> g_taxa;

    aw::LCAmapping s2i_map;     //supertree to input tree mapping
    aw::LCAmapping i2s_map;     //input tree to supertree mapping
    
    aw::TreetaxaMap s_nmap;        //taxamap for super tree    
    std::vector<aw::LCA> g_lca;   //lca's for input trees
    aw::LCA s_lca;   //lca's for supertree      

    {   //reading inputs
        const std::string filename = trees_filename;
        { // read trees
            std::ifstream ifs;
            std::istream &is = filename.empty() ? std::cin : ifs;
            if (!filename.empty()) {
                ifs.open(filename.c_str());
                if (!ifs) ERROR_exit("cannot read file '" << filename << "'");
            }            

            //first tree is supertree...
            MSG("Reading supertree: ");
            if (!aw::stream2tree(is, s_tree, s_taxa)) ERROR_exit("No supertree tree found in file '" << filename << "'");
            
            //reading other input trees
            MSG("Reading input trees: ");            
            for (;;) {                
                aw::Tree t;
                aw::idx2name t_names;                
                if (!aw::stream2tree(is, t, t_names)) break;
                g_taxa.push_back(t_names);
                g_trees.push_back(t);                
            }
        }
        MSG("Input trees: " << g_trees.size());
        if (g_trees.empty()) ERROR_exit("No input trees found in file '" << filename << "'");
    }

    // map input trees taxa labels
    aw::TaxaMap taxamap;  //to store all taxon and global id
    std::vector<aw::TreetaxaMap> g_nmaps; //for mapping taxamap and idx2name (of a tree)
    {
        g_nmaps.resize(g_taxa.size());
        for (unsigned int i=0,iEE=g_taxa.size(); i<iEE; ++i) {
            aw::idx2name &n = g_taxa[i];
            taxamap.insert(n);
            g_nmaps[i].create(n,taxamap);
        }
        MSG("Taxa: " << taxamap.size());
    }


    {   //map supertree labels
        aw::idx2name &n = s_taxa;       
        taxamap.insert(n);
        s_nmap.create(n,taxamap);
    }

    //strip extra leaves from the input trees for RF (opt - 3) and assessment (opt - 4)
    if(option==3 || option==4) {        
        for (unsigned int i=0,iEE=g_trees.size(); i<iEE; ++i) {
            int leaves0 = 0;
            int trim = 0;
            std::vector<unsigned int> leaves_trim;
            TREE_FOREACHLEAF(v,g_trees[i]) {
                leaves0++;
                if(g_trees[i].degree(v)!=1)  
                    MSG("ERROR");
                unsigned int gid_v = g_nmaps[i].gid(v);
                
                if(!s_nmap.exists(gid_v)) {                    
                    leaves_trim.push_back(v);
                }
            }
            if(leaves_trim.size()>0) 
                MSG("Tree: "<<i<<", trimming "<<leaves_trim.size());

            //actually trimming leaves
            BOOST_FOREACH(const unsigned int &c, leaves_trim) {                
                leaves0--;  trim++;
                if(g_trees[i].degree(c)!=1) {
                    MSG("Singleton node!");
                    continue;
                }
                if(!g_trees[i].trim_leaf(c))
                    ERROR_exit("Error trimming leaf...");
            }
        }        
    }

    if(option==1) {  //STATS OF ALL INPUT TREES .....
        MSG("Computing basic statistics...");
        unsigned int s_cnt=0;
        unsigned int unr_cnt=0;
        unsigned int r_cnt = 0;
        unsigned int no_phy_cnt = 0;
        unsigned int knee = 0;
        bool flag1 = true;
        unsigned int totIntNodes = 0;

        unsigned int avg_chld = 0;

        for (aw::Tree::iterator_dfs vl=s_tree.begin_dfs(),vEE=s_tree.end_dfs(); vl!=vEE; ++vl) {

             if(vl.idx==0) {
                 if(flag1) {
                     MSG("OpenTree root's degree: "<<s_tree.degree(vl.idx));
                     avg_chld += s_tree.degree(vl.idx);
                     flag1 = false;
                     totIntNodes++;
                 }
                 
                 continue;
             }
             if(vl.parent==0)
                 MSG(" Label of one of the three big clades:"<<s_taxa[vl.idx]);

             switch (vl.direction) {
                case aw::PREORDER: {
                    if(s_tree.is_leaf(vl.idx)) {
                        if(s_tree.degree(vl.parent)==2)
                            knee++;
                        s_cnt++;                     

                        continue;
                    }

                    if(s_tree.degree(vl.idx)>2) {
                        avg_chld += s_tree.degree(vl.idx) - 1;
                        totIntNodes++;
                    }
                    if(s_tree.degree(vl.idx)>3)
                        unr_cnt++;
                    else if(s_tree.degree(vl.idx)==2)
                        no_phy_cnt++;
                    else
                        r_cnt++;
                }
             }
        }
        output<<"OpenTree leaves: "<<s_cnt<<", unresolved nodes: "<<unr_cnt<<", resolved nodes: "<<r_cnt<<", degree two nodes: "<<no_phy_cnt<<", knees: "<<knee<<", avg_chld: "<<avg_chld/totIntNodes<<std::endl;

        std::vector<unsigned int> del_node;

        float totAvgChld = 0;
        for(unsigned int i=0; i<g_trees.size(); ++i) {
            unsigned int t_cnt=0;  //count leaves
            unsigned int unres_cnt=0;
            unsigned int r_cnt=0;
            unsigned int n_phy_cnt = 0;
            unsigned int knee = 0;
            
            bool flag1 = true;
            float totIntNodes = 0;
            unsigned int avg_chld = 0;
            
            for (aw::Tree::iterator_dfs vl=g_trees[i].begin_dfs(),vEE=g_trees[i].end_dfs(); vl!=vEE; ++vl) {
                 if(vl.idx==0) {
                     if(flag1) {
                        avg_chld += g_trees[i].degree(vl.idx);
                        flag1 = false;
                        totIntNodes++;
                     }
                     continue;
                 }

                 switch (vl.direction) {
                    case aw::PREORDER: {
                        if(g_trees[i].is_leaf(vl.idx)) {
                            if(g_trees[i].degree(vl.parent)==2 && vl.parent != g_trees[i].root)
                                knee++;
                            t_cnt++;
                            continue;
                        }
                        if(g_trees[i].degree(vl.idx)>2) {
                            avg_chld += g_trees[i].degree(vl.idx) - 1;
                            totIntNodes++;
                        }
                        if(g_trees[i].degree(vl.idx)>3) {
                            unres_cnt++;
                        } else if(g_trees[i].degree(vl.idx)==2) {
                            //del_node.push_back(vl.idx);
                            n_phy_cnt++;
                        }
                        else
                            r_cnt++;
                    }
                 }
            }

            BOOST_FOREACH(const unsigned int &j,del_node) {
                std::cout<<"^"<<j;
                g_trees[i].nphy_delete(j);
            }

            totAvgChld += avg_chld/totIntNodes;

            if(taxonomy && i==0)
                output<<" Taxonomy leaves: "<<t_cnt<<", unresolved nodes (non-root): "<<unres_cnt<<", resolved nodes  (non-root): "<<r_cnt<<", degree two nodes (non-root): "<<n_phy_cnt<<", knee: "<<knee<<", avg_chld: "<<avg_chld/totIntNodes<<std::endl;
            else
                output<<" Study tree "<<i<<" leaves: "<<t_cnt<<", unresolved nodes (non-root): "<<unres_cnt<<", resolved nodes  (non-root): "<<r_cnt<<", degree two nodes (non-root): "<<n_phy_cnt<<", knee: "<<knee<<", avg_chld: "<<avg_chld/totIntNodes<<std::endl;
        }

        output<<" Average node child: "<<totAvgChld/g_trees.size();

    }
    if(option==2) { //Counting Resolved Triplets ---------------------------------
        MSG("Counting resolved triplets...");
        long long int s_trip = 0;
        {   // compute cluster for the s_tree again
            TREE_POSTORDER2(v,s_tree) {
                if (s_tree.is_leaf(v.idx)){
                   s_tree.update_clst(v.idx,1);
                }
                else {
                    unsigned int count = 0;
                    BOOST_FOREACH(const unsigned int &c,s_tree.children(v.idx,v.parent)) {
                        count = count + s_tree.return_clstSz(c);
                    }
                    s_tree.update_clst(v.idx,count);                    
                }
            }
        }

        {   // compute cluster for the s_tree again
            TREE_POSTORDER2(v,s_tree) {
                if(s_tree.return_clstSz(v.idx) > 2) {
                    std::vector<unsigned int> ch_clst;
                    BOOST_FOREACH(const unsigned int &c,s_tree.children(v.idx,v.parent)) {
                        ch_clst.push_back(s_tree.return_clstSz(c));
                    }
                    for(int i=0; i<ch_clst.size(); ++i) {
                        if(ch_clst[i]<2) continue;
                        unsigned int sm = 0;
                        for(int j=0; j<ch_clst.size(); ++j) {
                            if(j==i) continue;
                            sm += ch_clst[j];
                        }
                        s_trip +=  (util::nchoosek(ch_clst[i],2))*sm;                        
                    }
                }
            }
        }
        output<<"\nTriplets in OpenTree: "<<s_trip;   

        //triplets of all input trees
        for(int m=0; m<g_trees.size(); ++m) {
            long long int s_trip = 0;
            {   // compute cluster for the s_tree again
                TREE_POSTORDER2(v,g_trees[m]) {
                    if (g_trees[m].is_leaf(v.idx)){
                       g_trees[m].update_clst(v.idx,1);
                    }
                    else {
                        unsigned int count = 0;
                        BOOST_FOREACH(const unsigned int &c,g_trees[m].children(v.idx,v.parent)) {
                            count = count + g_trees[m].return_clstSz(c);
                        }
                        g_trees[m].update_clst(v.idx,count);
                    }
                }
            }

            {   // compute cluster for the s_tree again
                TREE_POSTORDER2(v,g_trees[m]) {
                    if(g_trees[m].return_clstSz(v.idx) > 2) {
                        std::vector<unsigned int> ch_clst;
                        BOOST_FOREACH(const unsigned int &c,g_trees[m].children(v.idx,v.parent)) {
                            ch_clst.push_back(g_trees[m].return_clstSz(c));
                        }
                        for(int i=0; i<ch_clst.size(); ++i) {
                            if(ch_clst[i]<2) continue;
                            unsigned int sm = 0;
                            for(int j=0; j<ch_clst.size(); ++j) {
                                if(j==i) continue;
                                sm += ch_clst[j];
                            }
                            s_trip +=  (util::nchoosek(ch_clst[i],2))*sm;
                        }
                    }
                }
            }
            output<<"\nTriplets in Tree "<<m<<": "<<s_trip;
       }
    } 

    else if(option == 3) { //Computing RF distance
        MSG("Computing RF distance...");
       //for g_trees count clades ....    
        std::vector<unsigned int> g_nodes;
        {
            for(unsigned int i=0; i<g_trees.size(); ++i) {
                unsigned int cst = 0;
                TREE_POSTORDER2(v,g_trees[i]) {
                    if(g_trees[i].root!=v.idx && g_trees[i].degree(v.idx)>2) //counting unique and non-trivial clades
                        cst++;                   
                }
                g_nodes.push_back(cst);
            }
        } 

       {  //store LCAs
            aw::LCA lca;
            for (unsigned int i=0,iEE=g_trees.size(); i<iEE; ++i) {
                lca.create(g_trees[i]);
                g_lca.push_back(lca); }
            s_lca.create(s_tree);
        }

        aw::SubtreeParent<aw::Tree>  s_parent;
        s_parent.create(s_tree);
        bool tax_case = false;  //flag for taxonomy or studies
        
        unsigned int avg_ATE = 0;

        //main loop---------------------------------------------------
        for (unsigned int i=0,iEE=g_trees.size(); i<iEE; ++i) {
            MSG("Tree:"<<i);
            if(!taxonomy)
                output<<"Studies Tree: "<<i<<std::endl;

            int g_leaves = 0;
            TREE_FOREACHLEAF(v,g_trees[i]) {                
                    g_leaves++;
            }

            if(taxonomy)
                output<<" Taxonomy_leaves: "<<g_leaves;
            else
                output<<" Stree_leaves: "<<g_leaves;

            //set flag for taxonomy/studies
            if(taxonomy && i==0)  tax_case = true;
            else tax_case = false;

            std::queue<unsigned int> s_nodes; //nodes of s_tree in post-order
            unsigned int rooti;

            aw::SubtreeParent<aw::Tree>  g_parent;
            g_parent.create(g_trees[i]);

            //LCA for leaf nodes
            i2s_map.update_LCA_leaves(s_nmap,g_nmaps[i],s_tree,g_trees[i],tax_case);
            if(tax_case)
                s2i_map.update_LCA_leaves(g_nmaps[i],s_nmap,g_trees[i],s_tree,tax_case);
            else
                s2i_map.update_LCA_leaves_rel(g_nmaps[i],s_nmap,g_trees[i],s_tree,i2s_map, s_parent);

            unsigned int s_nds = 0; //counting internal nodes in the s_tree relative to a tree
            {   //Computing cluster size for supertree: computed based on leaf mapping
                unsigned int count;

                unsigned int irrelevent = 0;
                TREE_POSTORDER2(v,s_tree) {
                    if (!s_tree.is_leaf(v.idx)) {   //internal node case
                        if(!tax_case) {   //study trees case
                            if((s2i_map.mapping(v.idx) != NONODE) && (s2i_map.mapping(v.idx) != NONODE - 1))  //g_leaf mapped to it
                                s_tree.update_clst(v.idx,1);
                            else {
                                count = 0;
                                unsigned int pp = 0;
                                BOOST_FOREACH(const unsigned int &c,s_tree.children(v.idx,v.parent)) {
                                    count = count + s_tree.return_clstSz(c);
                                    if(s_tree.return_clstSz(c)>0) pp++;
                                }
                                s_tree.update_clst(v.idx,count);                                
                                if(count > 1 && count!=g_leaves && pp>1)  {
                                    s_nodes.push(v.idx);  //add only non-trivial clusters
                               }
                                if(count > 1 && pp>1) {
                                    rooti = v.idx;
                                }
                            }
                        }
                        else {   //taxonomy case
                            count = 0;
                            unsigned int pp = 0;
                            BOOST_FOREACH(const unsigned int &c,s_tree.children(v.idx,v.parent)) {
                                count = count + s_tree.return_clstSz(c);
                                if(s_tree.return_clstSz(c)>0) pp++;
                            }
                            s_tree.update_clst(v.idx,count);
                            if(count > 1  && s_tree.root!=v.idx && pp>1) s_nodes.push(v.idx);
                            else irrelevent++;
                        }
                    }
                    else {   //laef node case
                        if(s2i_map.mapping(v.idx)!=NONODE) {
                            s_tree.update_clst(v.idx,1);                            
                        }
                        else {
                            if(tax_case) ERROR_exit("Taxonomy must have all leaves!");
                            s_tree.update_clst(v.idx,0);
                        }
                    }
                }

                s_nds = s_nodes.size(); // all relevent nodes                
            }

            {   // compute cluster & score for input tree nodes....
                TREE_POSTORDER2(v,g_trees[i]) {
                    g_trees[i].update_score(v.idx,0);
                    if (g_trees[i].is_leaf(v.idx)){
                        unsigned int ggid = g_nmaps[i].gid(v.idx);
                        if(s_nmap.exists(ggid)) g_trees[i].update_clst(v.idx,1);
                        else g_trees[i].update_clst(v.idx,0);
                    }
                    else {
                        unsigned int count = 0;
                        BOOST_FOREACH(const unsigned int &c,g_trees[i].children(v.idx,v.parent)) {
                            count = count + g_trees[i].return_clstSz(c);
                        }

                        g_trees[i].update_clst(v.idx,count);
                    }
                }
            }


            //update LCA mapping for internal nodes
            s2i_map.update_LCA_internals(g_lca[i],s_tree,s_nmap,tax_case,'s');
            i2s_map.update_LCA_internals(s_lca,g_trees[i],g_nmaps[i],tax_case,'g');

            //unsigned int ints_nodes = s_nodes.size() - 1; //substracting root node
            output<<" OpenTree_Clades_Assessed: "<<s_nodes.size();

            unsigned int spt = 0;

            // for all s_tree nodes in the queue....................................
            while (!s_nodes.empty()) {
                unsigned int s_node = s_nodes.front();  s_nodes.pop();

                if(s_tree.root == s_node)  continue;
                unsigned int g_map = s2i_map.mapping(s_node);

                //FOR RF----
                if(s_tree.return_clstSz(s_node)== g_trees[i].return_clstSz(g_map))
                    spt++;
                
            }

            if(spt>g_nodes[i]) ERROR_exit("spt:"<<spt<<" g_nodes[i]:"<<g_nodes[i]);            
            output<<" Supported_Clades: "<<spt;
            output<<" STree_Int_Nodes: "<<g_nodes[i];
            output<<" RF_Distance: "<<(s_nds - spt + g_nodes[i]-spt)<<std::endl;
            if(s_nds+g_nodes[i] > 0)
                output<<" ATE: "<< ((s_nds - spt + g_nodes[i]-spt)*100)/(s_nds+g_nodes[i]);
            else
                output<<" ATE: "<<0;
            output<<std::endl;
            if(s_nds+g_nodes[i] > 0)
                avg_ATE += ((s_nds - spt + g_nodes[i]-spt)*100)/(s_nds+g_nodes[i]);
            else
                avg_ATE += 0;
            output<<"-----------------------------------------------------"<<std::endl;
        }
        output<<"Average dissimilarity: "<<avg_ATE/g_trees.size();
    }
    else if(option == 4) { //Assessing trees...
        MSG("Assessing trees following wilkinson et al.")

        {   //store lcas
            aw::LCA lca;
            for (unsigned int i=0,iEE=g_trees.size(); i<iEE; ++i) {
                lca.create(g_trees[i]);
                g_lca.push_back(lca); }
            s_lca.create(s_tree);
        }

        unsigned int Novel_nd1, Novel_nd2, Novel_nd3;
        {
            // compute cluster for the s_tree again
            TREE_POSTORDER2(v,s_tree) {
                if (s_tree.is_leaf(v.idx)){
                   s_tree.update_clst(v.idx,1);
                }
                else {
                    unsigned int count = 0;
                    BOOST_FOREACH(const unsigned int &c,s_tree.children(v.idx,v.parent)) {
                        count = count + s_tree.return_clstSz(c);
                    }
                    s_tree.update_clst(v.idx,count);
                }
            }
        }

        MSG("After LCA!");

        int * support;
        support = new int [s_tree.node_size()];
        for(int i=0,iEE=s_tree.node_size(); i<iEE; ++i)
            support[i] = 0;

        float * ws;
        ws = new float [s_tree.node_size()];
        for(int i=0,iEE=s_tree.node_size(); i<iEE; ++i)
            ws[i] = 0.0;

        int * ss;
        ss = new int [s_tree.node_size()];
        for(int i=0,iEE=s_tree.node_size(); i<iEE; ++i)
            ss[i] = 0;

        int * conflict;
        conflict = new int [s_tree.node_size()];
        for(int i=0,iEE=s_tree.node_size(); i<iEE; ++i)
            conflict[i] = 0;

        int * permit;
        permit = new int [s_tree.node_size()];
        for(int i=0,iEE=s_tree.node_size(); i<iEE; ++i)
            permit[i] = 0;

        aw::SubtreeParent<aw::Tree>  s_parent;
        s_parent.create(s_tree);

        bool tax_case = false;  //flag for taxonomy or studies

        //main loop---------------------------------------------------
        for (unsigned int i=0,iEE=g_trees.size(); i<iEE; ++i) {
            MSG("Tree:"<<i);

            int g_leaves = 0;
            TREE_FOREACHLEAF(v,g_trees[i]) {
                if(g_trees[i].degree(v)==1)
                    g_leaves++;
            }

            //set flag for taxonomy/studies
            if(taxonomy && i==0)  tax_case = true;
            else tax_case = false;

            std::queue<unsigned int> s_nodes; //nodes of s_tree in post-order

            aw::SubtreeParent<aw::Tree>  g_parent;
            g_parent.create(g_trees[i]);

            //LCA for leaf nodes

            i2s_map.update_LCA_leaves(s_nmap,g_nmaps[i],s_tree,g_trees[i],tax_case);
            if(tax_case)
                s2i_map.update_LCA_leaves(g_nmaps[i],s_nmap,g_trees[i],s_tree,tax_case);
            else
                s2i_map.update_LCA_leaves_rel(g_nmaps[i],s_nmap,g_trees[i],s_tree,i2s_map, s_parent);

            {   //Computing cluster size for supertree: computed based on leaf mapping
                unsigned int count;
                unsigned int irrelevent = 0;
                TREE_POSTORDER2(v,s_tree) {
                    if (!s_tree.is_leaf(v.idx)) {   //internal node case
                        if(!tax_case) {   //study trees case
                            if((s2i_map.mapping(v.idx) != NONODE) && (s2i_map.mapping(v.idx) != NONODE - 1))  //g_leaf mapped to it
                                s_tree.update_clst(v.idx,1);
                            else {
                                count = 0;
                                BOOST_FOREACH(const unsigned int &c,s_tree.children(v.idx,v.parent))
                                    count = count + s_tree.return_clstSz(c);
                                s_tree.update_clst(v.idx,count);
                                if(count > 1 && s_tree.root != v.idx)  s_nodes.push(v.idx);  //add only non-trivial clusters
                            }
                        }
                        else {   //taxonomy case
                            count = 0;
                            BOOST_FOREACH(const unsigned int &c,s_tree.children(v.idx,v.parent))
                                count = count + s_tree.return_clstSz(c);
                            s_tree.update_clst(v.idx,count);
                            if(count > 1  && s_tree.root != v.idx) s_nodes.push(v.idx);
                            else irrelevent++;
                        }
                    }
                    else {   //laef node case
                        if(s2i_map.mapping(v.idx)!=NONODE) {
                            s_tree.update_clst(v.idx,1);
                        }
                        else {
                            if(tax_case) ERROR_exit("Taxonomy must have all leaves!");
                            s_tree.update_clst(v.idx,0);
                        }
                    }
                }
                MSG("\nIrrelevent: "<<irrelevent);
            }

            {   // compute cluster & score for input tree nodes....
                TREE_POSTORDER2(v,g_trees[i]) {
                    g_trees[i].update_score(v.idx,0);
                    if (g_trees[i].is_leaf(v.idx)){
                        unsigned int ggid = g_nmaps[i].gid(v.idx);
                        if(s_nmap.exists(ggid)) g_trees[i].update_clst(v.idx,1);
                        else g_trees[i].update_clst(v.idx,0);
                    }
                    else {
                        unsigned int count = 0;
                        BOOST_FOREACH(const unsigned int &c,g_trees[i].children(v.idx,v.parent)) {
                            count = count + g_trees[i].return_clstSz(c);
                        }
                        g_trees[i].update_clst(v.idx,count);
                    }
                }
            }


            //update LCA mapping for internal nodes
            s2i_map.update_LCA_internals(g_lca[i],s_tree,s_nmap,tax_case,'s');
            i2s_map.update_LCA_internals(s_lca,g_trees[i],g_nmaps[i],tax_case,'g');

            //unsigned int ints_nodes = s_nodes.size() - 1; //substracting root node
            MSG("Total Clades Assessed: "<<s_nodes.size());

            // for all s_tree nodes in the queue....................................
            while (!s_nodes.empty()) {
                unsigned int s_node = s_nodes.front();  s_nodes.pop();

                if(s_tree.root == s_node)  continue;

                unsigned int g_map = s2i_map.mapping(s_node);

                //MSG(s_tree.return_clstSz(s_node)<<"-"<<g_trees[i].return_clstSz(g_map));

                if(s_tree.return_clstSz(s_node) < 2 || s_tree.return_clstSz(s_node) == g_leaves) {   // IRRELEVENT ....
                    //MSG(node<<" irrelevent node");
                }
                else if(s_tree.return_clstSz(s_node)== g_trees[i].return_clstSz(g_map)) {    // SUPPORT .....
                    support[s_node]++;
                    int b;
                    if(g_trees[i].return_score(g_map)==0) {
                        b = 1;
                        bool flag = true;
                        unsigned int node_par = s_parent.parent(s_node);
                        while(flag) {
                            if(s2i_map.mapping(node_par)==g_map)
                                b++;
                            else {
                                flag = false;
                                continue;
                            }
                            node_par = s_parent.parent(node_par);
                        }
                        g_trees[i].update_score(g_map,b);  // score stores how many nodes of s_tree lca to it
                    }
                    else {
                        b = g_trees[i].return_score(g_map);
                    }

                    ws[s_node] += 1.0/b;
                    if(b==1)  ss[s_node]++;
                }
                else {   // CONFLICT or PERMIT ......
                    std::vector<unsigned int> gmap_child;
                    unsigned int gmap_par = g_parent.parent(g_map);
                    g_trees[i].children(g_map,gmap_par,gmap_child);

                    if(gmap_child.size()>2) {
                        int cnt = 0;

                        BOOST_FOREACH(const unsigned int &c,g_trees[i].children(g_map,gmap_par)) {
                            unsigned int cmap = i2s_map.mapping(c);
                            unsigned int lcapar = s_lca.lca(cmap,s_node);

                            if(lcapar==s_node)
                                cnt = cnt + g_trees[i].return_clstSz(c);
                        }

                        if(cnt==s_tree.return_clstSz(s_node)) {
                            permit[s_node]++;     //MSG(node<<" consistent")
                        } else {
                            conflict[s_node]++;    //MSG(node<<" conflict");
                        }
                    }
                    else {
                        conflict[s_node]++;      //MSG(node<<" conflict");
                    }
                }
            }
        }

        float avg_V = 0.0;
        float avg_Vplus = 0.0;
        float avg_Vminus = 0.0;
        float cnt = 0.0;
        float avg_ws = 0.0;
        float avg_ss = 0.0;

        output <<"<--- ASSESSMENT OF SUPERTREE CLADES --->"<< std::endl<<std::endl;

        output <<"V"<<"\t"<<"V+"<<"\t"<<"V-"<<"\t"<<"ws"<<"\t"<<"ss"<< std::endl;

        {   // compute cluster for the s_tree again
            TREE_POSTORDER2(v,s_tree) {
                if (s_tree.is_leaf(v.idx)){
                   s_tree.update_clst(v.idx,1);
                }
                else {
                    unsigned int count = 0;
                    BOOST_FOREACH(const unsigned int &c,s_tree.children(v.idx,v.parent)) {
                        count = count + s_tree.return_clstSz(c);
                    }
                    s_tree.update_clst(v.idx,count);
                }
            }
        }

        //summarizing results now...
        unsigned int support_count = 0;
        TREE_POSTORDER2(v,s_tree) {
            if((s_tree.return_clstSz(v.idx)!=1) && (s_tree.root!=v.idx)) {                

                cnt += 1.0;
                if(support[v.idx]>=2)
                    support_count++;
                //output <<support[v.idx]<<"\t"<<conflict[v.idx]<<"\t"<<permit[v.idx]<<"\t"<<ws[v.idx]<<"\t"<<ss[v.idx]<< std::endl;
                float V = 0.0;
                float V_plus = 0.0;
                float V_minus = 0.0;
                if((support[v.idx]+conflict[v.idx])!=0)
                    V = ((support[v.idx] - conflict[v.idx])*1.0)/(support[v.idx] + conflict[v.idx]);
                if((support[v.idx]+conflict[v.idx]+permit[v.idx])!=0) {
                    V_plus = ((support[v.idx] - conflict[v.idx] + permit[v.idx])*1.0)/(support[v.idx] + conflict[v.idx] + permit[v.idx]);
                    V_minus = ((support[v.idx] - conflict[v.idx] - permit[v.idx])*1.0)/(support[v.idx] + conflict[v.idx] + permit[v.idx]);
                }
                avg_V += V;
                avg_Vplus += V_plus;
                avg_Vminus += V_minus;
                avg_ws += ws[v.idx];
                avg_ss += ss[v.idx];
                output <<std::setprecision(2) << std::fixed <<V<<"\t"<<std::setprecision(2) << std::fixed <<V_plus<<"\t"<<std::setprecision(2) << std::fixed <<V_minus<<"\t"<<std::setprecision(2) << std::fixed <<ws[v.idx]<<"\t"<<std::setprecision(2) << std::fixed <<ss[v.idx]<< std::endl;
                
            }
        }

        MSG("Nodes that are supported by two or more trees: "<<support_count);

        std::cout<<"\nClades assessed: "<<cnt;
        output<<"\nTree Assessment: average across clades --->"<< std::endl<< std::endl;
        output <<"Avg_V"<<"\t"<<"Avg_V+"<<"\t"<<"Avg_V-"<<"\t"<<"Avg_ws"<<"\t"<<"Avg_ss"<< std::endl;
        output<<std::setprecision(4) << std::fixed <<avg_V/cnt<<"\t"<<std::setprecision(4) << std::fixed <<avg_Vplus/cnt<<"\t"<<std::setprecision(4) << std::fixed <<avg_Vminus/cnt<<"\t"<<std::setprecision(4) << std::fixed <<avg_ws/cnt<<"\t"<<std::setprecision(4) << std::fixed <<avg_ss/cnt;

        delete [] support;
        delete [] conflict;
        delete [] permit;
        delete [] ws;
        delete [] ss;

        MSG("\nDone assessing!!!")
    }


    return (EXIT_SUCCESS);
}


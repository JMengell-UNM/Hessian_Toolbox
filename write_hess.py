def Wriste_Hess(self, opp_dict=None):
    
    temp_dict={}
    for key in self.hess_dict.keys():
        if opp_dict!=None and key in opp_dict.keys():
            temp_dict[key]=opp_dict[key]
        else:
            temp_dict[key]=self.hess_dict[key]
            
    print(temp_dict)
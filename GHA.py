import numpy
import random

from itertools import accumulate
import itertools


from tqdm import tqdm

count_all=0
class GA:

    def __init__(self,functionname="Himmelblau"):
        self.bounds = {"Himmelblau": {'x': [-5, 5], 'y': [-5, 5]},
                       "McCormick": {'x': [-1.5, 4], 'y': [-3, 3]},
                       "Bohachevsky": {'x': [-10, 10], 'y': [-10, 10]},
                       "otherfunction": {'x1': [0, 0], 'x2': [0, 0], 'x3': [0, 0]},
                       }
        self.functionname=functionname
        self.population_bin=[]
        self.population_bin_sep=[]
        self.population_val=[]
        self.l_sep=[]
        self.l_sum=0
        self.bound=self.bounds[functionname]
        self.father=[]
        self.son=[]
        self.popsize=100
        self.precision=0.01
        self.pc = 0.8  # 交叉概率
        self.pm = 0.2  # 变异概率
        self.z_list=[]
        self.z_f = []
        self.z_p_sum=[]

        self.epoch = 100
    def Himmelblau(self, x, y):
        return (x ** 2 + y - 11) ** 2 + (x + y ** 2 - 7) ** 2
    def McCormick(self, x, y):
        return numpy.sin(x + y) + (x - y) ** 2 - 1.5 * x + 2.5 * y + 1
    def Bohachevsky(self, x, y):
        return x ** 2 + 2 * y ** 2 - 0.3 * numpy.cos(3 * numpy.pi * x) * numpy.cos(4 * numpy.pi * y) + 0.3

    def CalculateBinLenth(self):
        print("===========================计算个体基因长度=============================")
        bound=self.bounds[self.functionname]
        for key, value in bound.items():
            print('变量: ', key, '区间: ', value)
            w=(value[1]-value[0])/self.precision
            m=0
            while not (2 ** (m - 1) < w<= 2 ** m):
                m = m + 1
            self.l_sep.append(m)
            print('精度为 {0} 时，变量 {1} 在区间 {2} 需要编码的数的个数为 {3} 个'.format(self.precision,key,value,w))
            print("其二进制编码长度为: {0}".format(m))
            print("--------------------------------------------------------------")
            self.l_sum=self.l_sum+m
        print("编码总长度为{0}".format(self.l_sum))

    def PopRandomInit(self):
        print("===========================种群随机初始化===============================")
        print("种群中的个体：")
        for i in range(self.popsize):
            individual = ''.join([random.choice('01') for j in range(self.l_sum)])
            self.population_bin.append(individual)
        count = 0  # 设置初始计数
        for k in range(len(self.population_bin)):
            print(self.population_bin[k], end=' ')
            count += 1  # 开始计数
            if count % 5 == 0:  # 每5个换行
                print(end='\n')
        print("--------------------------------------------------------------")
        print("个体数量一共: {0} 个".format(self.popsize))
        print("--------------------------------------------------------------")
        print("个体基因长度为: {0}".format(self.l_sum))
        print("--------------------------------------------------------------")

    def Slection(self):

        for i in range(len(self.population_bin)):
            code_sep1=self.population_bin[i][0:self.l_sep[0]]
            code_sep2=self.population_bin[i][self.l_sep[0]:]
            temple=[]
            temple.append(code_sep1)
            temple.append(code_sep2)
            self.population_bin_sep.append(temple)
            temple=[]
        print("种群中的个体（x,y分开写）：")
        print(self.population_bin_sep)
        print("--------------------------------------------------------------")


        for i in range(len(self.population_bin_sep)):
            temple = []
            temple.append(int(self.population_bin_sep[i][0], 2) / (2 ** len(self.population_bin_sep[i][0]) - 1) * (self.bound['x'][1]- self.bound['x'][0]) +self.bound['x'][0])
            temple.append(int(self.population_bin_sep[i][1], 2) / (2 ** len(self.population_bin_sep[i][1]) - 1) * (self.bound['y'][1]- self.bound['y'][0]) +self.bound['y'][0])
            self.population_val.append(temple)
        print("二进制串对应的十进制的值为 ： ")
        print(self.population_val)
        print("--------------------------------------------------------------")


        z_p=[]
        if self.functionname=="Himmelblau":
            for i in self.population_val:
                z= self.Himmelblau(i[0],i[1])
                self.z_list.append(z)
                self.z_f.append(1/z*1000)
            print("种群个体十进制值带入 Himmelblau 函数对应的函数值: \n{0}".format(self.z_list))
            print("--------------------------------------------------------------")
            print("种群个体适应度适应度为:（乘1000取倒数）\n {0}".format(self.z_f))
            print("--------------------------------------------------------------")
        if self.functionname=="McCormick":
            for i in self.population_val:
                z= self.McCormick(i[0],i[1])
                self.z_list.append(z)
                self.z_f.append(1/z*1000)
            print("种群个体十进制值带入 McCormick 函数对应的函数值: \n{0}".format(self.z_list))
            print("--------------------------------------------------------------")
            print("种群个体适应度适应度为:（乘1000取倒数）\n {0}".format(self.z_f))
            print("--------------------------------------------------------------")

        for i in self.z_f:
            z_f_sum=sum(self.z_f)
            temple=i/z_f_sum
            z_p.append(temple)
        print("种群适应度概率为: \n{0}".format(z_p))
        print("--------------------------------------------------------------")
        self.z_p_sum = list(accumulate(z_p))
        print("种群适应度累积概率和为: \n{0}".format(self.z_p_sum))


        print("================================选择模块================================")
        count = 0
        count_temp=1
        father1 = "0"
        father2 = "1"
        while count < self.popsize // 2 and father1 != father2:
            print("第{0}次选择".format(count_temp))
            r1 = random.random()
            print("选择概率r1= {0}".format(r1))
            for i in range(len(self.z_p_sum)):
                if r1 <= self.z_p_sum[i]:
                    father1 = self.population_bin[i]
                    print("选择第{0}个个体：{1}".format(i+1,father1))
                    break
                if self.z_p_sum[i] < r1 <= self.z_p_sum[i + 1]:
                    father1 = self.population_bin[i + 1]
                    print("选择第{0}个个体：{1}".format(i+2, father1))
                    break

            r2 = random.random()
            print("选择概率r2= {0}".format(r2))
            for i in range(len(self.z_p_sum)):
                if r2 <= self.z_p_sum[i]:
                    father2 = self.population_bin[i]
                    print("选择第{0}个个体：{1}".format(i+1, father2))
                    break
                if self.z_p_sum[i] < r2 <= self.z_p_sum[i + 1]:
                    father2 = self.population_bin[i + 1]
                    print("选择第{0}个个体：{1}".format(i+2, father2))
                    break
            if father1 != father2:
                self.father.append([father1, father2])
                count += 1
            else:
                print("选择的两个个体重复！！！重新选择！！！！")
                father1 = "0"
                father2 = "1"
            count_temp+=1
            print("--------------------------------------------------------------")
        print("选择完毕!!!!!!\n共选择 {0} 对父本，分别是\n{1}".format(count,self.father))


        self.population_bin.clear()
        self.population_bin_sep.clear()
        self.population_val.clear()
        self.z_list.clear()
        self.z_f.clear()
        self.z_p_sum.clear()


    def Crossover(self):
        print("================================开始交叉================================")
        for i in range(len(self.father)):

            r=random.randint(1, self.l_sum-1)
            print("生成随机断点{0}".format(r))

            rc=random.random()
            print("生成随机概率rc:",rc)
            if rc<self.pc:
                son1=self.father[i][0][:r]+self.father[i][1][r:]
                son2 = self.father[i][1][:r] + self.father[i][0][r:]
                self.son.append([son1,son2])
                print("rc<pc,进行交叉")
                print("--------------------------------------------------------------")
            else:
                self.son.append(self.father[i])
                print("rc>pc,不进行交叉")
                print("--------------------------------------------------------------")
        print("交叉后的子代:\n{0}".format(self.son))
        self.father=self.son.copy()
        self.son.clear()
    def Mutation(self):
        print("================================开始变异================================")
        count=1
        global count_all
        self.father=list(itertools.chain(*self.father))
        print("父代：",self.father)

        for i in range(len(self.father)):
            r=random.randint(1, self.l_sum)
            rm = random.random()
            print("第{0}次变异，生成随机变异点{1},生成随机概率:{2}".format(count,r,rm),end=" ")
            if rm<self.pm:
                temp = self.father[i][:r-1] + ('1' if self.father[i][r-1] == '0' else '0') + self.father[i][r:]
                self.son.append(temp)
                print("rm<pm,进行变异")
                print("--------------------------------------------------------------")
            else:
                self.son.append(self.father[i])
                print("rm>pm,不进行变异")
                print("--------------------------------------------------------------")
            count+=1
        print("变异后的子代:\n{0}".format(self.son))
        self.population_bin = self.son.copy()
        self.father.clear()
        self.son.clear()


        for i in range(len(self.population_bin)):
            code_sep1=self.population_bin[i][0:self.l_sep[0]]
            code_sep2=self.population_bin[i][self.l_sep[0]:]
            temple=[]
            temple.append(code_sep1)
            temple.append(code_sep2)
            self.population_bin_sep.append(temple)
            temple=[]
        print("种群中的个体（x,y分开写）：")
        print(self.population_bin_sep)
        print("--------------------------------------------------------------")


        for i in range(len(self.population_bin_sep)):
            temple = []
            temple.append(int(self.population_bin_sep[i][0], 2) / (2 ** len(self.population_bin_sep[i][0]) - 1) * (self.bound['x'][1]- self.bound['x'][0]) +self.bound['x'][0])
            temple.append(int(self.population_bin_sep[i][1], 2) / (2 ** len(self.population_bin_sep[i][1]) - 1) * (self.bound['y'][1]- self.bound['y'][0]) +self.bound['y'][0])
            self.population_val.append(temple)
        print("二进制串对应的十进制的值为 ： ")
        print(self.population_val)
        print("--------------------------------------------------------------")


        z_p=[]
        if self.functionname=="Himmelblau":
            for i in self.population_val:
                z= self.Himmelblau(i[0],i[1])
                self.z_list.append(z)
            print("种群个体十进制值带入 Himmelblau 函数对应的函数值: \n{0}".format(self.z_list))
            print("--------------------------------------------------------------")
        count_all+=1
        if count_all == self.epoch :
            best=min(self.z_list)
            print("{0} 轮迭代完成后，最优值为{1}".format(self.epoch,best))
        self.population_bin_sep.clear()
        self.population_val.clear()
        self.z_list.clear()



ga1=GA()
ga1.CalculateBinLenth()
ga1.PopRandomInit()
for i in tqdm(range(ga1.epoch)):
    print("================================================================")
    print("================================================================")
    print("第{0}次迭代".format(i+1))
    print("================================================================")
    print("================================================================")
    ga1.Slection()
    ga1.Crossover()
    ga1.Mutation()

# -*- coding: utf-8 -*-
"""
Created on Sun Oct 28 23:22:10 2018

@author: yangcy
"""

import quandl
import os
import math
import copy
import datetime
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

class Implementation(object):
    def __init__(self, cycle, group_name):
        self.cycle = cycle
        self.save_path = "C:/Users/yangcy/Desktop/2018 fall/hedge funds/pre1/data/" + group_name + "_" + str(cycle) + "/"
        self.gether = pd.Series(name = group_name, index = ["IR","rtn","std","max drawdown"])
    
    def downloadData(self):
#        #CBOT Soybean
#        self.s1 = quandl.get("CHRIS/CME_S1", authtoken="jtRZ1ZYz35_a28X4dJ6U", start_date="2002-01-28", end_date="2018-10-26")
#        self.s10 = quandl.get("CHRIS/CME_S10", authtoken="jtRZ1ZYz35_a28X4dJ6U", start_date="2002-01-28", end_date="2018-10-26")
#        #Soybean Oil
#        self.bo1 = quandl.get("CHRIS/CME_BO1", authtoken="jtRZ1ZYz35_a28X4dJ6U", start_date="2002-01-28", end_date="2018-10-26")
#        self.bo10 = quandl.get("CHRIS/CME_BO10", authtoken="jtRZ1ZYz35_a28X4dJ6U", start_date="2002-01-28", end_date="2018-10-26")
#        
        #Live Cattle
        self.lc1 = quandl.get("CHRIS/CME_lc1", authtoken="jtRZ1ZYz35_a28X4dJ6U", start_date="2002-01-28", end_date="2018-10-26")
        self.lc7 = quandl.get("CHRIS/CME_lc7", authtoken="jtRZ1ZYz35_a28X4dJ6U", start_date="2002-01-28", end_date="2018-10-26")
        #Feeder Cattle
        self.fc1 = quandl.get("CHRIS/CME_fc1", authtoken="jtRZ1ZYz35_a28X4dJ6U", start_date="2002-01-28", end_date="2018-10-26")
        self.fc7 = quandl.get("CHRIS/CME_fc7", authtoken="jtRZ1ZYz35_a28X4dJ6U", start_date="2002-01-28", end_date="2018-10-26")
        
        #Brent Crude
        self.ice_b1 = quandl.get("CHRIS/ICE_b1", authtoken="jtRZ1ZYz35_a28X4dJ6U", start_date="2002-01-28", end_date="2018-10-26")
        self.ice_b12 = quandl.get("CHRIS/ICE_b13", authtoken="jtRZ1ZYz35_a28X4dJ6U", start_date="2002-01-28", end_date="2018-10-26")
        #NYMEX WTI Crude
        self.cme_cl1 = quandl.get("CHRIS/CME_cl1", authtoken="jtRZ1ZYz35_a28X4dJ6U", start_date="2002-01-28", end_date="2018-10-26")
        self.cme_cl12 = quandl.get("CHRIS/CME_cl13", authtoken="jtRZ1ZYz35_a28X4dJ6U", start_date="2002-01-28", end_date="2018-10-26")
        
#        #NYMEX Gasoline
#        self.cme_rb1 = quandl.get("CHRIS/CME_rb1", authtoken="jtRZ1ZYz35_a28X4dJ6U", start_date="2002-01-28", end_date="2018-10-26")
#        self.cme_rb12 = quandl.get("CHRIS/CME_rb12", authtoken="jtRZ1ZYz35_a28X4dJ6U", start_date="2002-01-28", end_date="2018-10-26")
#        #ICE Gas Oil
#        self.ice_g1 = quandl.get("CHRIS/ICE_g1", authtoken="jtRZ1ZYz35_a28X4dJ6U", start_date="2002-01-28", end_date="2018-10-26")
#        self.ice_g12 = quandl.get("CHRIS/ICE_g12", authtoken="jtRZ1ZYz35_a28X4dJ6U", start_date="2002-01-28", end_date="2018-10-26")
#        #NYMEX Heating Oil
#        self.cme_ho1 = quandl.get("CHRIS/CME_ho1", authtoken="jtRZ1ZYz35_a28X4dJ6U", start_date="2002-01-28", end_date="2018-10-26")
#        self.cme_ho12 = quandl.get("CHRIS/CME_ho12", authtoken="jtRZ1ZYz35_a28X4dJ6U", start_date="2002-01-28", end_date="2018-10-26")
#        
        return
    
    def mergeDate(self, data_list):
        data_list_new = list()
        inter = list()
        for i in range(len(data_list)):
            index_list = list(data_list[i].index)
            if i == 0:
                inter = index_list
            else:
                inter = list(set(inter) & set(index_list))
        inter.sort()
        for i in range(len(data_list)):
            data_list_new.append(data_list[i].loc[inter])
        return data_list_new
    
    def calYt(self, settle_price_all):
        Yt = pd.DataFrame(columns = list(range(0,len(settle_price_all.columns),2)))
        for i in range(0,len(settle_price_all.columns),2):
            Yt[i] = (settle_price_all[i])/(settle_price_all[i+1])
        return Yt
        
    def calWeight(self, Yt):
        n = len(Yt.columns)
        Yt_rank = Yt.rank(axis = 1)
        weight = pd.DataFrame(0, index = Yt.index, columns = Yt.columns)
        weight[Yt_rank == n] = 1
        weight[Yt_rank == 1] = -1
        return weight
    
    def calNetValue(self, group1):
        settle_price_all = pd.DataFrame(columns = list(range(len(group1))))
        for i in range(len(group1)):
            settle_price_all[i] = group1[i]['Settle']
        settle_price_all = settle_price_all.fillna(method = 'ffill')
        
        Yt = self.calYt(settle_price_all)
        weight = self.calWeight(Yt)
        
        drop_list = list(range(1,len(group1),2))
        settle_price = settle_price_all.drop(drop_list, axis = 1)
        
        self.fwdrtn = settle_price.shift(-1)/settle_price - 1
        
        df = pd.DataFrame(index = settle_price.index)
        df['cumcount'] = range(len(df))
        df['period'] = np.floor(df.cumcount / self.cycle)
        df['day'] = df.cumcount - df.period * self.cycle
        self.full_date_index = df.index
        self.trade_date_index = df.loc[df.day == 0].index
        
        weight_real = pd.DataFrame(np.nan, index = settle_price.index, columns = settle_price.columns)
        weight_real.loc[self.trade_date_index] = weight
        if self.cycle != 1:
            self.weight = weight_real.ffill(limit = self.cycle - 1)#如果cycle为1就会出问题
        else:
            self.weight = weight_real
        
        portfolio_rtn = (1. + self.fwdrtn).groupby(df['period']).cumprod()
        portfolio_rtn = self.weight * portfolio_rtn
        portfolio_rtn = portfolio_rtn.sum(axis = 1).shift(1).fillna(0.)
        portfolio_rtn = portfolio_rtn + 1

        temp = portfolio_rtn.loc[self.trade_date_index]
        temp = temp.cumprod()

        net_value =  pd.Series(np.nan, index = self.full_date_index)
        net_value[self.trade_date_index] = temp
        if self.cycle != 1:
            net_value = net_value.ffill(limit = self.cycle - 1) #如果cycle为1就会出问题

        portfolio_rtn.loc[self.trade_date_index] = 1.
        self.net_value_without_fee = net_value * portfolio_rtn
        
        return
    
    def drawMaxDrawdown(self):
        drawdown = []
        for i in range(len(self.full_date_index)):
            di = self.lsdata.iloc[i]
            dwi = (di/(self.lsdata[:(i+1)]).max())-1
            drawdown.append(dwi)
        drawdown = pd.Series(drawdown,index = self.full_date_index)*100
        date1 = pd.to_datetime(self.full_date_index,format = '%Y%m%d')
        plt.plot(date1,drawdown)
        #IR, MEAN, STDY
        plt.ylabel('drawdown(%)')
        plt.title('drawdown')
        plt.savefig(os.path.join(self.save_path, 'drawdown'))   
        #plt.show()
        plt.close()
        return (drawdown).min()
    
    def drawWealthCurve(self, win, year):
        win = 240
        self.net_rtn_without_fee = self.net_value_without_fee / self.net_value_without_fee.shift(1) - 1#
        self.net_rtn_without_fee = self.net_rtn_without_fee.fillna(0.)
        lsnet = copy.deepcopy(self.net_rtn_without_fee)
        lsdata = (lsnet+1).cumprod()
        self.lsdata = lsdata
        self.lsnet = lsnet
        date1 = pd.to_datetime(lsdata.index,format = '%Y%m%d')
        drw = self.drawMaxDrawdown()
        
        plt.plot(date1,lsdata)
        #IR, MEAN, STDY
        plt.ylabel('Wealth Curve')
        plt.title('Wealth Curve')
        length = lsdata.max()-lsdata.min()
        
        # 计算一共有多少天，除去开头和结尾的 NA
        temp = pd.DataFrame()
        temp['lsnet'] = lsnet
        temp['count'] = 1
        temp['count'] = temp['count'].cumsum()
        temp = temp.dropna()
        temp = temp.iloc[-1, 1] - temp.iloc[0, 1] # column 1 refers to ['count']
        
        
        text1 = "IR: "+ str(round(((self.net_rtn_without_fee).mean()/np.nanstd(self.net_rtn_without_fee))*math.sqrt(win),2))
        text2 = "rtn: " + str(round((100*(lsdata.iloc[-1])**(win*1.0/temp)-100),2)) +" %"
        text3 = "std: " + str(round((np.std(100*lsnet)*math.sqrt(win)),2)) +" %"
#        text4 = "turnover: " + str(round(self.Turnover.mean(),2))
        text5 = "max drawdown: " + str(round(drw,2)) + " %"

        self.gether['IR'] = text1[4:]
        self.gether['rtn'] = text2[5:]
        self.gether['std'] = text3[5:]
#        self.gether['turnover'] = text4[10:]
        self.gether['max drawdown'] = text5[14:]
        
#        plt.text(datetime.datetime(year, 4, 4, 0, 0), (lsdata).max(), text1)  
#        plt.text(datetime.datetime(year, 4, 4, 0, 0), (lsdata).max()-0.11*length, text2)  
#        plt.text(datetime.datetime(year, 4, 4, 0, 0), (lsdata).max()-0.22*length, text3)
#        plt.text(datetime.datetime(year, 4, 4, 0, 0), (lsdata).max()-0.24*length, text4)
#        plt.text(datetime.datetime(year, 4, 4, 0, 0), (lsdata).max()-0.33*length, text5)
        
        plt.savefig(os.path.join(self.save_path, 'Wealth Curve'))   
        #plt.show()
        plt.close()
        return self.gether
    
    
if __name__ == '__main__':
    cycle = 5
    backtest = Implementation(cycle, 'Crude Oils New')
    print("Downloading Data")
    backtest.downloadData()
    
    group1 = list([backtest.ice_b1, backtest.ice_b12, backtest.cme_cl1, backtest.cme_cl12])
#    group1 = list([backtest.lc1, backtest.lc7, backtest.fc1, backtest.fc7])
#    group1 = list([backtest.s1, backtest.s10, backtest.bo1, backtest.bo10])
#    group1 = list([backtest.cme_rb1, backtest.cme_rb12, backtest.ice_g1, backtest.ice_g12, backtest.cme_ho1, backtest.cme_ho12])
    
    print("Calculating")
    group1 = backtest.mergeDate(group1)
    backtest.calNetValue(group1)
    result = backtest.net_value_without_fee
    
    print("Plotting")
    gether = backtest.drawWealthCurve(240, 2002)
    gether.to_csv(os.path.join(backtest.save_path, 'result.csv'))



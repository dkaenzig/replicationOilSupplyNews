%% Load the data

if strcmp(dataFrequency,'M')
    [dataRaw,labels] = xlsread('data/DataUS.xlsx','Export');    % load the data set
    header           = labels(1,2:end);
    datesStringRaw   = labels(2:end,1);
    datesNumRaw      = (str2double(datesStringRaw{1}(1:4))+(str2double(datesStringRaw{1}(end-1:end))-1)*1/12: ...
                        1/12:str2double(datesStringRaw{end}(1:4))+(str2double(datesStringRaw{end}(end-1:end))-1)*1/12)';

    Tall = size(dataRaw,1);

    % generate variables for series contained in data
    for ii = 1:length(header)
        eval([header{ii} '= dataRaw(:,ii);'])
    end

    % Overview of series %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % FRED: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % WTISPLC: Spot Crude Oil Price: West Texas Intermediate (WTI)
    % POILWTIUSDM: Global price of WTI Crude
    % DCOILWTICO: Crude Oil Prices: West Texas Intermediate (WTI) - Cushing, Oklahoma
    % POILBREUSDM: Global price of Brent Crude
    % DCOILBRENTEU: Crude Oil Prices: Brent - Europe
    % BRENTEXT: Crude Oil Prices: Brent - Europe (extended using WTI)
    % WTIBRENT: WTI price(extended by Brent from 2005M1 because it is a better proxy for the global oil price)
    % INDPRO: Industrial Production Index
    % IPMAN: Industrial Production: Manufacturing (NAICS)
    % AWHMAN: Average Weekly Hours of Production and Nonsupervisory Employees: Manufacturing
    % UNRATE: Civilian Unemployment Rate
    % PAYEMS: All Employees: Total Nonfarm Payrolls
    % CPIAUCSL: Consumer Price Index for All Urban Consumers: All Items
    % CPILFESL: Consumer Price Index for All Urban Consumers: All Items Less Food and Energy
    % CPIENGSL: Consumer Price Index for All Urban Consumers: Energy
    % CUSR0000SAN: Consumer Price Index for All Urban Consumers: Nondurables
    % CUSR0000SAD: Consumer Price Index for All Urban Consumers: Durables	
    % CUSR0000SAS: Consumer Price Index for All Urban Consumers: Services
    % PCEPI: Personal Consumption Expenditures: Chain-type Price Index
    % PCEPILFE: Personal Consumption Expenditures Excluding Food and Energy (Chain-Type Price Index)
    % DNRGRG3M086SBEA: Personal consumption expenditures: Energy goods and services (chain-type price index)
    % DDURRG3M086SBEA: Personal consumption expenditures: Durable goods (chain-type price index)
    % DNDGRG3M086SBEA: Personal consumption expenditures: Nondurable goods (chain-type price index)
    % DSERRG3M086SBEA: Personal consumption expenditures: Services (chain-type price index)
    % PPIACO: Producer Price Index for All Commodities
    % WPSFD49207: Producer Price Index by Commodity for Final Demand: Finished Goods
    % WPSFD4131: Producer Price Index by Commodity for Final Demand: Finished Goods Less Foods and Energy 
    % WPU0561: Producer Price Index by Commodity for Fuels and Related Products and Power: Crude Petroleum (Domestic Production)
    % A576RC1: Compensation of Employees, Received: Wage and Salary Disbursements
    % A132RC1: Compensation of Employees, Received: Wage and Salary Disbursements: Private Industries
    % FF: Effective Federal Funds Rate
    % DGS1:	1-Year Treasury Constant Maturity Rate
    % DGS10: 10-Year Treasury Constant Maturity Rate
    % IPG211111CS: Industrial Production: Mining: Crude oil
    % MICH: University of Michigan: Inflation Expectation
    % PCE: Personal Consumption Expenditures
    % PCEDG: Personal Consumption Expenditures: Durable Goods
    % PCEND: Personal Consumption Expenditures: Nondurable Goods
    % PCES: Personal Consumption Expenditures: Services
    % DNRGRC1M027SBEA: Personal consumption expenditures: Energy goods and services
    % VXOCLS: CBOE S&P 100 Volatility Index: VXO
    % VXOCLSEXT: CBOE S&P 100 Volatility Index: VXO (extended using Bloom's series based on actual volatilities)
    % VIXCLS: CBOE Volatility Index: VIX
    % NNUSBIS: Narrow Effective Exchange Rate for United States
    % RNUSBIS: Real Narrow Effective Exchange Rate for United States
    % NBUSBIS: Broad Effective Exchange Rate for United States
    % RBUSBIS: Real Broad Effective Exchange Rate for United States
    % TWEXMMTH: Trade Weighted U.S. Dollar Index: Major Currencies
    % TWEXMPA: Real Trade Weighted U.S. Dollar Index: Major Currencies
    % TWEXBMTH: Trade Weighted U.S. Dollar Index: Broad
    % TWEXBPA: Real Trade Weighted U.S. Dollar Index: Broad
    % TWEXOMTH:	Trade Weighted U.S. Dollar Index: Other Important Trading Partners
    % TWEXOPA: Real Trade Weighted U.S. Dollar Index: Other Important Trading Partners
    % Datastream: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % EIA1955: Crude Oil Prod. World (Thousand Barrels / Day) - TURNOVER BY VOLUME
    % EIA1944: Crude Oil Total Production OPEC - TURNOVER BY VOLUME	
    % EIA1954: Crude Oil Total Prod Non-OPEC - TURNOVER BY VOLUME	
    % EIA2004: Crude Oil Production United States - TURNOVER BY VOLUME
    % WDPQPTLMP: WD PETROLEUM: PROD,OIL SUPPLY-CRUDE OIL INCL LEASE CONDENSATE
    % USINCOXRP: US STOCKS - CRUDE OIL EXCL STRATEGIC PETROLEUM RESERVES VOLN
    % EIA1533: Crude Oil Stocks Total (million barrels) - INVENTORY VOLUME US
    % EIA1532: Crude Oil Stocks Non-SPR MMBBL - INVENTORY VOLUME US
    % EIA1531: Crude Oil Stocks Strategic Petroleum Reserve (SPR) - INVENTORY VOLUME US (Million Barrels)
    % WDSCFIIWP: WD CRUDE OIL & LIQUID FUELS: INVENTORY NET WITHDRAWALS,STOCK DRA	
    % OCSCFICP:	OECD CRUDE OIL & LIQUID FUELS: INVENTORIES-COMMERCIAL INVENTORY
    % EIA1976: Petroleum Total Stocks OECD (million barrels) - INVENTORY VOLUME (includes both crude oil and petroleum products (refined crude oil))
    % EIA1976EXT: Petroleum Total Stocks OECD - INVENTORY VOLUME (extended as in Kilian and Murphy)
    % EIA1541: Petroleum Total Stocks US (million barrels) - INVENTORY VOLUME
    % SPCOMP: S&P 500 COMPOSITE - PRICE INDEX
    % SPINDS: S&P INDUSTRIAL - PRICE INDEX
    % SP5EIND: S&P500 ES INDUSTRIALS - PRICE INDEX
    % SP5EUTL: S&P500 ES UTILITIES - PRICE INDEX
    % SP5GTRS: S&P500 TRANSPORTATION - PRICE INDEX
    % SP5EENE: S&P500 ES ENERGY - PRICE INDEX
    % DJINDUS: DOW JONES INDUSTRIALS - PRICE INDEX
    % DJCMP65: DOW JONES COMPOSITE 65 STOCK AVE - PRICE INDEX
    % DJUTILS: DOW JONES UTIILITIES - PRICE INDEX
    % DJTRSPT: DOW JONES TRANSPORTATION - PRICE INDEX
    % NASCOMP: NASDAQ COMPOSITE - PRICE INDEX
    % NYSEALL: NYSE COMPOSITE - PRICE INDEX 
    % (all stock prices also in EOP format, name: XXXXEOP)
    % USCOCODMA: US REFINERS ACQUISITION COST OF DOMESTIC CRUDE OIL CURN	
    % USCOCOIMA: US REFINERS ACQUISITION COST OF IMPORTED CRUDE OIL CURN	
    % USCOCOA: US REFINERS ACQUISITION COST OF DOM. & IMPORTED CRUDE OIL CURN
    % CRUDOILEXT: CRUDE OIL SPOT PRICE (EOP starting from 1983, before
    % average)
    % OILGSUS: US-DS Oil & Gas - PRICE INDEX	
    % OILGPUS: US-DS Oil & Gas Prod - PRICE INDEX	
    % OILESUS: US-DS Oil/Eq Svs/Dst - PRICE INDEX
    % BMATRUS: US-DS Basic Materials - PRICE INDEX	
    % BRESRUS: US-DS Basic Resource - PRICE INDEX
    % INDMTUS: US-DS Ind. Met & Mines - PRICE INDEX	
	% MNINGUS: US-DS Mining - PRICE INDEX
    % INDUSUS: US-DS Industrials - PRICE INDEX	
    % INDGSUS: US-DS Industrial Goods & Services - PRICE INDEX
    % AERSPUS: US-DS Aero/Defence - PRICE INDEX	
    % INDTRUS: US-DS Inds Transpt - PRICE INDEX
    % CNSMGUS: US-DS Consumer Gds - PRICE INDEX
    % AUTMBUS: US-DS Auto & Parts - PRICE INDEX
    % FDBEVUS: US-DS Food & Bev - PRICE INDEX
    % PERHHUS: US-DS Pers & H/H Gds - PRICE INDEX
    % HLTHCUS: US-DS Health Care - PRICE INDEX	
    % CNSMSUS: US-DS Consumer Svs - PRICE INDEX	
    % RTAILUS: US-DS Retail - PRICE INDEX
    % TRLESUS: US-DS Travel & Leis - PRICE INDEX
    % TELCMUS: US-DS Telecom - PRICE INDEX
    % UTILSUS: US-DS Utilities - PRICE INDEX
    % FINANUS: US-DS Financials - PRICE INDEX	
    % INSURUS: US-DS Insurance - PRICE INDEX	
    % RLESTUS: US-DS Real Estate - PRICE INDEX	
    % FINSVUS: US-DS Financial Svs(3) - PRICE INDEX	
    % TECNOUS: US-DS Technology - PRICE INDEX
    % TECHDUS: US-DS Tch H/W & Eq - PRICE INDEX
    % TOTMKUS: US-DS Market - PRICE INDEX								
	% AEROSUS: US-DS Aerospace - PRICE INDEX
    % AIRLNUS: US-DS Airlines - PRICE INDEX
    % AUTOSUS: US-DS Automobiles - PRICE INDEX
    % BIOTCUS: US-DS Biotechnology - PRICE INDEX
    % CHMCLUS: US-DS Chemicals - PRICE INDEX
    % CHEMSUS: US-DS Commodity Chem - PRICE INDEX
    % DEFENUS: US-DS Defense - PRICE INDEX
    % DELSVUS: US-DS Delivery Svs - PRICE INDEX
    % ELECTUS: US-DS Electricity - PRICE INDEX
    % INDTRUS: US-DS Inds Transpt - PRICE INDEX
    % MARINUS: US-DS Marine Transpt - PRICE INDEX
    % TRUCKUS: US-DS Trucking - PRICE INDEX
    % TRNSVUS: US-DS Transpt Svs - PRICE INDEX
    % PHRMCUS: US-DS Pharm - PRICE INDEX
    % RAILSUS: US-DS Railroads - PRICE INDEX
    % ARCAOIL: NYMEX oil index
    % UXOM:	EXXON MOBIL - PRICE INDEX			
    % RDSB:	ROYAL DUTCH SHELL B - PRICE INDEX
    % RDSA:	ROYAL DUTCH SHELL A - PRICE INDEX
    % BP: BP - PRICE INDEX
    % USTOTPRCF: US TERMS OF TRADE REBASED TO 1975=100 NADJ	US 
    % USEXPTOTB: TOTAL EXPORTS ON A BALANCE OF PAYMENTS BASIS CURA	
    % USIMPTOTB: US TOTAL IMPORTS ON A BALANCE OF PAYMENTS BASIS CURA		
	% USBALTOTB: US GOODS & SERVICES BALANCE ON A BALANCE OF PAYMENTS BASIS CURA			
    % USEXPGDSB: US EXPORTS F.A.S. CURA	
    % USIMPGDSB: US IMPORTS F.A.S. CURA	
    % USVISGDSB: US VISIBLE TRADE BALANCE F.A.S.-F.A.S. CURA
    % CRSP: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % MDRET: Value-Weighted Return-incl. dividends	
    % MRET:	Value-Weighted Return-excl. dividends
    % SPRET:	Return on the S&P 500 Index
    % SP: Level of the S&P 500 Index
    % Misc: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % GLOBACTOLD: Kilians old global activity measure
    % GLOBACTNEW: Kilians corrected global activity measure
    % OECDIP: Industrial production of OECD countries
    % EBPG: Gilchrists excess bond premium
    % GZSPREAD: Gilchrist's GZ spread
    % OECDINVCH: OECD oil inventory change from Kilian and Murphy
    % OECDSTOCKSKilian: OECD above ground oil inventories as constructed in Kilian and Murphy
    % OECDSTOCKSKilianSA: OECD above ground oil inventories as constructed in Kilian and Murphy (SA)
    % OECD6MNEIP: OECD+6NME Industrial Production from Baumeister Hamilton
    % WTIFUTXEOP: WTI futures price (X-month contract, EOP)
    % WTIFUTXAVG: WTI futures price (X-month contract, monthly average)
    % GPR: Caldara and Iacoviello geopolitical risk index (baseline)	
    % GPR_THREAT: Caldara and Iacoviello geopolitical risk index (threats)	
    % GPR_ACT: Caldara and Iacoviello geopolitical risk index (acts)	
    % GPR_BROAD: Caldara and Iacoviello geopolitical risk index (broad)	
    % GPR_NARROW: Caldara and Iacoviello geopolitical risk index (narrow)	
    % GPR_SAUDI_ARABIA: Caldara and Iacoviello geopolitical risk index SAUDI ARABIA
    % GPRH: Historical Caldara and Iacoviello geopolitical risk index (baseline)	
    % GPRHT: Historical Caldara and Iacoviello geopolitical risk index (threats)
    % GPRHA: Historical Caldara and Iacoviello geopolitical risk index (acts)
    % EPUNEWS: Economic policy uncertainty index (news based, Baker Bloom Davis)
    % BKEXPXM: Baumeister and Kilian oil price expectation XM (X=3 6 9 12)
    % BKEXP3MEXT, BKEXP6MEXT: Expectations extended with futures data
    % RRShock: Romer and Romer MP shock, extended by SilviaS
    % IFS (billateral exchange rates) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % USAUS	USBRA USCAN	USCHINA	USCOL USCZECH USDEN	USEGYPT USEUR USHUN USICELAND USINDIA USINDO USIRAN 
    % USIRAQ USISRAEL USJAP	USKOR USKUWAIT USLIBYA USMAL USMEX USNEWZ USNIG USNOR USOMAN USPHIL USPOL 
    % USQATAR USRUS USRUSF USSAUDI USSOUTHAF USSWE USSWI USTURK USEMIR USUK USVEN
    % IRAQ and KUWAIT have missing values
    % OECD IP: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % AUTIP	BELIP BRAIP CANIP CHLIP	COLIP CZEIP DEUIP DNKIP EA19IP ESPIP ESTIP EU28IP	
    % FINIP	FRAIP G7IP GBRIP GRCIP HUNIP INDIP IRLIP ISLIP ISRIP ITAIP JPNIP KORIP	
    % LTUIP	LUXIP LVAIP	MEXIP NLDIP	NORIP OECDIP OECDEIP POLIP PRTIP RUSIP SVKIP	
    % SVNIP	SWEIP TURIP	USAIP

    % define/construct additional series %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    OECDSTOCKS = OECDSTOCKSKilianSA; % (EIA1533./EIA1541.*EIA1976EXT;) alternatively use seasonally adjusted series
    GLOBACT = GLOBACTNEW; % use new (corrected) GLOBACT measure from Kilian
    
    % % compare my measure to Kilian and Murphy's
    % figure 
    % hold on
    % plot(OECDINV)
    % plot(OECDSTOCKS)
elseif strcmp(dataFrequency,'Q')
    [dataRaw,labels] = xlsread('data/DataUSQ.xlsx','Export');    % load the data set
    header           = labels(1,2:end);
    datesStringRaw   = labels(2:end,1);
    datesNumRaw      = (str2double(datesStringRaw{1}(1:4))+(str2double(datesStringRaw{1}(6))-1)*0.25: ...
                        0.25:str2double(datesStringRaw{end}(1:4))+(str2double(datesStringRaw{end}(6))-1)*0.25)';

    Tall = size(dataRaw,1);

    % generate variables for series contained in data
    for ii = 1:length(header)
        eval([header{ii} '= dataRaw(:,ii);'])
    end

    % Overview of series %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % FRED: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % WTISPLC: Spot Crude Oil Price: West Texas Intermediate (WTI)
    % GDPC1: Real Gross Domestic Product
    % CPIAUCSL: Consumer Price Index for All Urban Consumers: All Items
    % MICH: University of Michigan: Inflation Expectation
    % PCECC96: Real Personal Consumption Expenditures					
    % GPDIC1: Real Gross Private Domestic Investment
    % GCEC1: Real Government Consumption Expenditures and Gross Investment
    % NETEXP: Net Exports of Goods and Services
    % GDP: Gross Domestic Product
    % GDPDEF: Gross Domestic Product: Implicit Price Deflator
    % PRS85006022: Nonfarm Business Sector: Average Weekly Hours
    % Datastream: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % EIA1955: Crude Oil Prod. World MBBL/DAY - TURNOVER BY VOLUME
    % EIA1944: Crude Oil Total Production OPEC - TURNOVER BY VOLUME	
    % EIA1954: Crude Oil Total Prod Non-OPEC - TURNOVER BY VOLUME	
    % EIA2004: Crude Oil Production United States - TURNOVER BY VOLUME
    % USBALGDSB: US BOP BALANCE: GOODS CURA
    % USCURBALB: US CURRENT ACCOUNT BALANCE CURA
    % Misc: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % GLOBACTOLD: Kilians old global activity measure
    % GLOBACTNEW: Kilians corrected global activity measure
    % CPI6: SPF median inflation expectations (1 year horizon)
    % OECDSTOCKSEOP: monthly OECD oil stocks (as constructed above) aggregated by EOP
    % OECDSTOCKSAVG: monthly OECD oil stocks (as constructed above) aggregated by AVG
    % OECDSTOCKSSAEOP: monthly OECD oil stocks (as constructed above, SA) aggregated by EOP
    % OECDSTOCKSSAAVG: monthly OECD oil stocks (as constructed above, SA) aggregated by AVG
    % OECD6MNEIP: OECD+6NME Industrial Production from Baumeister Hamilton (aggregated by taking sum)
    
    % define/construct additional series %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    GLOBACT = GLOBACTNEW; % use new (corrected) GLOBACT measure from Kilian
end

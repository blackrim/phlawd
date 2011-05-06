#include "mainwindow.h"
#include "ui_mainwindow.h"
#include <stdlib.h>
#include <iostream>
using namespace std;

#include <QMainWindow>
#include <QFileDialog>
#include <QtSql>

MainWindow::MainWindow(QWidget *parent) :    QMainWindow(parent),ui(new Ui::MainWindow){
	ui->setupUi(this);
}

MainWindow::~MainWindow(){
	delete ui;
}

/*
  could allow for multiple connections but not doing that for now
 */
void MainWindow::load_sequencesqlite(){
	if (sdb.connectionName().length() != 0){
		cout << "existing connection" << endl;
		QSqlDatabase::database(sdb.databaseName()).close();
		QSqlDatabase::removeDatabase(sdb.databaseName());
	}
	QString filename = QFileDialog::getOpenFileName( this, tr("Open Document"), QDir::currentPath(), tr("Database files (*.db *.sql *.sqlite);;All files (*.*)"), 0, QFileDialog::DontUseNativeDialog );
	cout << filename.toStdString() << endl;
	sdb = QSqlDatabase::addDatabase("QSQLITE");
	sdb.setDatabaseName(filename);
	if (!sdb.open())
			cout << "open failed" << endl;
	else{
		ui->statusBar->showMessage("sequence database opened");
		ui->actionDB_stats->setEnabled(true);
		ui->actionDB_stats->setText("get stats: "+filename);
		ui->statusBar->showMessage("current: "+filename);
	}
	QSqlQuery q;
	cout << q.exec("SELECT Count(*) FROM sequence;") << endl;
	if (q.next()){
		QString numvalues (q.value(0).toString());
		cout << q.value(0).toString().toStdString() << endl;
		ui->textBrowser->setText("Number of sequences: "+numvalues);
	}
}

void MainWindow::get_sequencedb_stats(){
	QSqlQuery q;

	cout << q.exec("SELECT Count(*) FROM sequence;") << endl;
	if (q.next()){
		cout <<	 q.value(0).toString().toStdString()<< endl;
	}
	ui->textBrowser->setText("");
}

from flask_wtf import FlaskForm
from wtforms import StringField, PasswordField, FloatField, SubmitField
from wtforms.validators import DataRequired, Length

class MoleculeForm(FlaskForm):
    name = StringField("Name", validators=[DataRequired()])
    formula = StringField("Formula", validators=[DataRequired()])
    molecular_weight = FloatField("Molecular Weight", validators=[DataRequired()])
    structure = StringField("Structure (SMILES)", validators=[DataRequired()])
    submit = SubmitField("Add Molecule")

class LoginForm(FlaskForm):
    username = StringField("Username", validators=[DataRequired(), Length(min=3, max=20)])
    password = PasswordField("Password", validators=[DataRequired()])
    submit = SubmitField("Login")

class ReactionForm(FlaskForm):
    reaction = StringField('Reaction', validators=[DataRequired()])

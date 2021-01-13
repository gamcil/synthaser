import Select from 'react-select'

import { FilterList } from './filter/List'
import { RenameList } from './rename/List'

import { v4 as uuidv4 } from 'uuid'

export const RuleItem = props => {

  const handleUpdateRule = (label, value) => {
    props.handleChange({
      target: { name: label, value : value }
    })
  }

  const handleAddFilter = () => {
    handleUpdateRule(
      "filters",
      [{ uuid: uuidv4(), type: "", domains: [] }, ...props.data.filters]
    )
  }
  const handleChangeFilter = idx => events => {
    handleUpdateRule(
      "filters",
      props.data.filters.map((filter, _idx) => {
        if (_idx !== idx) return filter
        const data = { ...filter }
        events.forEach(event => {
          data[event.target.name] = event.target.value
        })
        return data
      })
    )
  }
  const handleRemoveFilter = idx => () => {
    handleUpdateRule(
      "filters",
      props.data.filters.filter((_, _idx) => _idx !== idx)
    )
  }

  const handleAddRename = () => {
    handleUpdateRule(
      "renames",
      [{ uuid: uuidv4(), "from": "", "after": [], "to": "" }, ...props.data.renames]
    )
  }
  const handleChangeRename = idx => event => {
    handleUpdateRule(
      "renames",
      props.data.renames.map((rename, _idx) => {
        if (_idx !== idx) return rename
        return { ...rename, [event.target.name]: event.target.value}
      })
    )
  }
  const handleRemoveRename = idx => () => {
    handleUpdateRule(
      "renames",
      props.data.renames.filter((_, _idx) => _idx !== idx)
    )
  }

  const domainOptions = props.domains.map(domain => ({
    label: domain.name,
    value: domain.name
  }))
  const handleChangeDomains = event => props.handleChange({
    target: {
      name: "domains",
      value: event ? event.map(e => e.value) : []
    }
  })

  const ruleOptions = props.rules
    .filter(rule => rule.uuid !== props.data.uuid)
    .map(rule => ({ label: rule.name, value: rule.uuid }))
  const handleChangeParent = event => props.handleChange({
    target: { name: "parent", value: event.value }
  })

  return (
    <li>
      <button
        type="button"
        onClick={props.handleRemove}
      >
        Delete
      </button>
      <div className="rule-field">
        <label htmlFor="ruleName">Name:</label>
        <input
          id="ruleName"
          type="text"
          name="name"
          value={props.data.name}
          onChange={props.handleChange}
        />
      </div>

      <div className="rule-field">
        <label htmlFor="ruleDomains">Domains:</label>
        <Select
          name="domains"
          options={domainOptions}
          onChange={handleChangeDomains}
          className="select"
          value={props.data.domains.map(d => domainOptions.find(o => o.label === d))}
          isMulti
        />
      </div>
      
      <div className="rule-field">
        <label htmlFor="ruleEvaluator">Evaluation expression:</label>
        <input
          id="ruleEvaluator"
          type="text"
          name="evaluator"
          value={props.data.evaluator}
          onChange={props.handleChange}
        />
      </div>

      <div className="rule-field">
        <label htmlFor="ruleParent">Parent rule:</label>
        <Select
          id="ruleParent"
          name="parent"
          className="select"
          options={ruleOptions}
          onChange={handleChangeParent}
          value={ruleOptions.find(o => o.label === props.data.parent)}
        />
      </div>

      <label htmlFor="ruleFilters">Domain filters:</label>
      <FilterList
        id="ruleFilters"
        rule={props.data}
        filters={props.data.filters}
        domains={props.domains}
        handleAdd={handleAddFilter}
        handleRemove={handleRemoveFilter}
        handleChange={handleChangeFilter}
      />

      <label htmlFor="ruleRename">Rename domains:</label>
      <RenameList
        id="ruleRename"
        rule={props.data}
        renames={props.data.renames}
        domains={props.domains}
        handleAdd={handleAddRename}
        handleRemove={handleRemoveRename}
        handleChange={handleChangeRename}
      />
    </li>
  )
}

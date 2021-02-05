import { useMemo } from 'react'
import React from 'react'
import Select from 'react-select'

import FilterList from './filter/List'
import RenameList from './rename/List'

import { v4 as uuidv4 } from 'uuid'


const RuleItem = ({
  name,
  domains,
  allDomains,
  evaluator,
  filters,
  renames,
  handleRemove,
  handleChange,
  handleUpdate,
  domainOptions,
}) => {

  const selectedDomains = useMemo(
    () => domains
      .map(d => ({
        label: d,
        innerValue: d,
        value: Math.random()
      })),
    [domains]
  )

  const handleChangeDomains = event => {
    handleUpdate("domains", event ? event.map(e => e.label) : [])
  }

  const handleAddFilter = () => {
    handleUpdate(
      "filters",
      [{ uuid: uuidv4(), type: "", domains: [] }, ...filters]
    )
  }
  const handleChangeFilter = uuid => events => {
    handleUpdate(
      "filters",
      filters.map(filter => {
        if (filter.uuid !== uuid) return filter
        const data = { ...filter }
        events.forEach(event => {
          data[event.target.name] = event.target.value
        })
        return data
      })
    )
  }
  const handleRemoveFilter = uuid => () => {
    handleUpdate(
      "filters",
      filters.filter(filter => filter.uuid !== uuid)
    )
  }

  const handleAddRename = () => {
    handleUpdate(
      "renames",
      [{
        uuid: uuidv4(),
        "from": "",
        "before": [],
        "after": [],
        "to": ""
      }, ...renames]
    )
  }
  const handleChangeRename = uuid => event => {
    handleUpdate(
      "renames",
      renames.map(rename => {
        if (rename.uuid !== uuid) return rename
        return { ...rename, [event.target.name]: event.target.value}
      })
    )
  }
  const handleRemoveRename = uuid => () => {
    handleUpdate(
      "renames",
      renames.filter(rename => rename.uuid !== uuid)
    )
  }

  return (
    <li>
      <button
        type="button"
        onClick={handleRemove}
      >
        Delete
      </button>

      <div className="rule-field">
        <label htmlFor="ruleName">Name:</label>
        <input
          id="ruleName"
          type="text"
          name="name"
          value={name}
          onChange={handleChange}
        />
      </div>

      <div className="rule-field">
        <label htmlFor="ruleDomains">Domains:</label>
        <Select
          name="domains"
          options={domainOptions}
          onChange={handleChangeDomains}
          className="select"
          value={selectedDomains}
          hideSelectedOptions={false}
          isMulti
        />
      </div>
      
      <div className="rule-field">
        <label htmlFor="ruleEvaluator">Evaluation expression:</label>
        <input
          id="ruleEvaluator"
          type="text"
          name="evaluator"
          value={evaluator}
          onChange={handleChange}
        />
      </div>

      <label htmlFor="ruleFilters">Domain filters:</label>
      <FilterList
        id="ruleFilters"
        filters={filters}
        domains={domains}
        allDomains={allDomains}
        handleAdd={handleAddFilter}
        handleRemove={handleRemoveFilter}
        handleChange={handleChangeFilter}
      />

      <label htmlFor="ruleRename">Rename domains:</label>
      <RenameList
        id="ruleRename"
        renames={renames}
        domains={domains}
        options={selectedDomains}
        handleAdd={handleAddRename}
        handleRemove={handleRemoveRename}
        handleChange={handleChangeRename}
        domainOptions={domainOptions}
      />
    </li>
  )
}

export default React.memo(RuleItem)
